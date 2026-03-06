#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.bam   = "./input/simplex.bam"
params.ref       = "./input/reference.fasta"
params.roi_bed   = "./input/UMIs_randomer_barcode.bed"
params.acc_r     = "./pairing_acc.R"
params.r10_data  = "./input/paired_output_R10.tsv"
params.hd_data   = "./input/paired_output_HD.tsv"
params.plot_r    = "./plot_acc.R"
params.mods_data = "./input/pairing_acc_TelN_moremods.tsv"

process bam2fastq {
    publishDir "data/"

    input:
    path bam

    output:
    path "simplex.fastq"

    script:
    """
    samtools fastq -T "*" $bam > simplex.fastq
    """
}

process minimap2_align {
    publishDir "data/"

    input:
    path fastq
    path ref

    output:
    path "minimap2_k5w3.bam"

    script:
    """
    minimap2 --MD -y -a $ref $fastq -k 5 -w 3 | \
      samtools sort -T tmp -o minimap2_k5w3.bam
    """
}

process filter_primary {
    publishDir "data/"

    input:
    path bam

    output:
    path "minimap2_pri.bam"

    script:
    """
    samtools view -bS -F 256 -F 16 -F 4 -F 2048 $bam | \
      samtools sort -T tMp -o minimap2_pri.bam ;
    """
}

process filter_concatemers {
    publishDir "data/"

    input:
    path bam

    output:
    path "noconca_RN.rn"
 
    script:
    """
    # The correct full size of a read including adapters should be 198 + 36 + 27 = 261, set threshold as 270 to allow small insertions
    samtools view $bam | \
      awk '{print \$1"\\t"\$4"\\t"length(\$10)}' | \
      awk '{print \$1"\\t"\$2+\$3}' | \
      awk '{if (\$2<=270) print \$0}' | \
      cut -f1 > noconca_RN.rn
    """
}

process self_annealing {
    publishDir "data/"

    input:
    path bam

    output:
    path "selfannealing_RN.rn"

    script:
    """
    samtools view $bam | \
      grep "GCGTGACATCGATGACACGC" | \
      cut -f1 > selfannealing_RN.rn
    """
}

process intersect_reads {
    publishDir "data/"

    input:
    path rn1
    path rn2

    output:
    path "non_conca_noSA_RN.rn"

    script:
    """
    cat $rn1 $rn2 $rn2 | sort | uniq -u > non_conca_noSA_RN.rn
    """
}

process filter_bam {
    publishDir "data/"

    input:
    path bam
    path rn

    output:
    path "filtered.bam"
    
    script:
    """
    samtools view -N $rn -o filtered.bam $bam
    """
}

process essential_regions {
    publishDir "data/"

    input:
    path bam

    output:
    path "filtered.bed"
    path "essential_RN.rn"

    script:
    """
    bedtools bamtobed -i $bam > filtered.bed
    
    # Filtered reads covering full barcodes on both strands
    less filtered.bed | \
      awk '{if (\$2 <= 44 && \$3 >= 226) print \$0}' | \
      cut -f4 > essential_RN.rn
    """
}

process final_filter_bam {
    publishDir "data/"

    input:
    path bam
    path bed
    path rn

    output:
    path "final_filterd.bam"
    
    script:
    """
    samtools view -N $rn -o final_filterd.bam $bam
    """
}

process index_filtered_bam {
    publishDir "data/"

    input:
    path bam

    output:
    path "final_filterd.bam.bai"

    script:
    """
    samtools index $bam
    """
}

process pileup {
    publishDir "data/"

    input:
    path bam
    path ref
    path roi_bed

    output:
    path "mpileup_output.maf"

    script:
    """
    samtools mpileup $bam -f $ref --no-output-ins --no-output-del --no-output-ends \
      -l $roi_bed -o mpileup_output.maf --output-QNAME -A -B --max-depth 0
    """
}

process split_bases {
    publishDir "data/base_information"

    input:
    path maf

    output:
    path "baseATpositions", type: 'dir'

    script:
    """
    mkdir -p baseATpositions

    awk '{print "echo "\$5">baseATpositions/"\$1"_"\$2"_a"}' $maf | sh
    awk '{print "echo "\$1"\\t"\$2"\\t"\$7">baseATpositions/"\$1"_"\$2"_b"}' $maf | sh
    """
}

process mutate_and_merge {
    publishDir "data/base_information"

    input:
    path baseATpos

    output:
    path "baseATpositions/merged_mutated", type: 'dir'

    script:
    """
    mkdir -p baseATpositions/a_mutated
    mkdir -p baseATpositions/b_mutated
    mkdir -p baseATpositions/merged_mutated

    for file in baseATpositions/*a; do
        name=\$(basename "\$file")
        awk -v FS="" '{for (i=1;i<=NF;i++) printf \$i"\\n"}' "\$file" |
        awk '{if (\$0 != "+" && \$0 != "-" && \$0 !~ /^[0-9]*\$/) print \$0}' \
        > baseATpositions/a_mutated/\${name}_mutated
    done

    for file in baseATpositions/*b; do
        name=\$(basename "\$file")
        awk '{print \$3}' "\$file" | sed 's/,/\\n/g' \
        > baseATpositions/b_mutated/\${name}_mutated
    done

    for file in baseATpositions/*a; do
        id=\$(basename "\$file" | awk -F'_' '{print \$1"_"\$2}')
        paste baseATpositions/a_mutated/\${id}_a_mutated \
            baseATpositions/b_mutated/\${id}_b_mutated \
        > baseATpositions/merged_mutated/\${id}_merged
    done
    """
}

process cal_pairing_acc {
    publishDir "data/", mode:"copy"

    input:
    path roi_bed
    path merged_dir
    path acc_r

    output:
    path "UMI_matching_TelN.tsv"

    script:
    """
    mkdir UMI_permut_bases_out 
    
    for r in \$(cut -f1 $roi_bed | sort -u); do
          Rscript $acc_r \$r ./
    done

    ls UMI_permut_bases_out/*raw*all.tsv | \
      awk -F'[/|_]' '{print \$0"\t"\$5}' | \
      while read a b;do awk '{print "'\$b'\t"\$0}' \$a; done | \
      grep -wv rn | grep -v "*" > UMI_matching_TelN.tsv
    """
}

process plot_acc {
    publishDir "data/", mode:"copy"

    input:
    path acc_teln
    path acc_r10
    path acc_hd
    path acc_mods
    path plot_acc_r

    output:
    path "UMI_matching.pdf"

    script:
    """
    grep -v "*" $acc_r10 > pairing_acc_R10.tsv
    grep -v "*" $acc_hd > pairing_acc_HD.tsv

    Rscript $plot_acc_r $acc_teln $acc_r10 $acc_hd $acc_mods UMI_matching.pdf 
    """
}

workflow {
    // Define channels for inputs
    bam_ch   = Channel.fromPath(params.bam)
    ref_ch   = Channel.fromPath(params.ref)   
    bed_ch   = Channel.fromPath(params.roi_bed)
    acc_r_ch = Channel.fromPath(params.acc_r)
    r10_ch   = Channel.fromPath(params.r10_data)
    hd_ch    = Channel.fromPath(params.hd_data)
    plot_r_ch = Channel.fromPath(params.plot_r)
    mods_ch  = Channel.fromPath(params.mods_data)

    // Stepwise pipeline
    fastq_ch     = bam2fastq(bam_ch)
    aln_ch       = minimap2_align(fastq_ch, ref_ch)
    primary_ch   = filter_primary(aln_ch)
    conca_ch     = filter_concatemers(primary_ch)
    selfann_ch   = self_annealing(primary_ch)
    inter_ch     = intersect_reads(conca_ch, selfann_ch)
    bam_filtered_ch      = filter_bam(primary_ch, inter_ch)
    essential_ch = essential_regions(bam_filtered_ch)
    bam_final_ch = final_filter_bam(bam_filtered_ch, essential_ch)
    index_ch     = index_filtered_bam(bam_final_ch)
    maf_ch       = pileup(bam_final_ch, ref_ch, bed_ch)
    bases_ch     = split_bases(maf_ch) 
    merged_ch    = mutate_and_merge(bases_ch)
    acc_ch       = cal_pairing_acc(bed_ch, merged_ch, acc_r_ch)
    plot_acc(acc_ch, r10_ch, hd_ch, mods_ch, plot_r_ch)
}
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.bam       = "./input/simplex.bam"
params.ref       = "./input/reference.fasta"
params.umi_bed   = "./input/UMI.bed"
params.ribo_bed  = "./input/ribo_trinu.bed"
params.umi_r_for = "./UMI_combining_bases_for.R"
params.umi_r_rev = "./UMI_combining_bases_rev.R"
params.umi_r_pairing = "./pairing_up_for_rev_UMI+channel.R"

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

process filter_primary_for {
    publishDir "data/"

    input:
    path bam

    output:
    path "minimap2_pri_forward.bam"

    script:
    """
    samtools view -bS -F 256 -F 16 -F 4 -F 2048 $bam | \
      samtools sort -T tMp -o minimap2_pri_forward.bam ;
    """
}

process filter_primary_rev {
    publishDir "data/"

    input:
    path bam

    output:
    path "minimap2_pri_reverse.bam"

    script:
    """
    samtools view -bS -F 256 -f 16 -F 4 -F 2048 $bam | \
      samtools sort -T tMp -o minimap2_pri_reverse.bam ;
    """
}

process filter_concatemers {
    publishDir "data/"

    input:
    path bam_for
    path bam_rev

    output:
    path "noconca_RN_for.rn"
    path "noconca_RN_rev.rn"
 
    script:
    """
    # The correct full size of a read including adapters should be 137 + 36 + 27 = 200, set threshold as 210 to allow small insertions
    samtools view $bam_for | \
      awk '{print \$1"\\t"\$4"\\t"length(\$10)}' | \
      awk '{print \$1"\\t"\$2+\$3}' | \
      awk '{if (\$2<=210) print \$0}' | \
      cut -f1 > noconca_RN_for.rn
    
    samtools view $bam_rev | \
      awk '{print \$1"\\t"\$4"\\t"length(\$10)}' | \
      awk '{print \$1"\\t"\$2+\$3}' | \
      awk '{if (\$2<=210) print \$0}' | \
      cut -f1 > noconca_RN_rev.rn
    """
}

process essential_regions {
    publishDir "data/"

    input:
    path bam_for
    path bam_rev

    output:
    path "essential_RN_for.rn"
    path "essential_RN_rev.rn"

    script:
    """
    bedtools bamtobed -i $bam_for | \
      awk '{if (\$2 <= 30 && \$3 >= 115) print \$0}' | \
      cut -f4 > essential_RN_for.rn

    bedtools bamtobed -i $bam_rev | \
      awk '{if (\$2 <= 30 && \$3 >= 115) print \$0}' | \
      cut -f4 > essential_RN_rev.rn
    """
}

process intersect_reads {
    publishDir "data/"

    input:
    path rn1_for
    path rn1_rev
    path rn2_for
    path rn2_rev

    output:
    path "non_conca_with_essential_RN_for.rn"
    path "non_conca_with_essential_RN_rev.rn"

    script:
    """
    cat $rn1_for $rn2_for | sort | uniq -d > non_conca_with_essential_RN_for.rn
    cat $rn1_rev $rn2_rev | sort | uniq -d > non_conca_with_essential_RN_rev.rn
    """
}

process filter_bam {
    publishDir "data/"

    input:
    path bam_for
    path bam_rev
    path rn_for
    path rn_rev

    output:
    path "filterd_for.bam"
    path "filterd_rev.bam"
    
    script:
    """
    samtools view -N $rn_for -o filterd_for.bam $bam_for
    samtools view -N $rn_rev -o filterd_rev.bam $bam_rev
    """
}

process index_filtered_bam {
    publishDir "data/"

    input:
    path bam_for
    path bam_rev

    output:
    path "filterd_for.bam.bai"
    path "filterd_rev.bam.bai"

    script:
    """
    samtools index $bam_for
    samtools index $bam_rev
    """
}

process pileup {
    publishDir "data/"

    input:
    path bam_for
    path bam_rev
    path ref
    path umi_bed

    output:
    path "mpileup_output_for_umi.maf"
    path "mpileup_output_rev_umi.maf"

    script:
    """
    samtools mpileup $bam_for -f $ref --no-output-ins --no-output-del --no-output-ends \
      -l $umi_bed -o mpileup_output_for_umi.maf --output-QNAME -A -B --max-depth 0

    samtools mpileup $bam_rev -f $ref --no-output-ins --no-output-del --no-output-ends \
      -l $umi_bed -o mpileup_output_rev_umi.maf --output-QNAME -A -B --max-depth 0
    """
}

process split_bases {
    publishDir "data/base_information"

    input:
    path maf_for
    path maf_rev

    output:
    path "baseATpositions_for", type: 'dir'
    path "baseATpositions_rev", type: 'dir'

    script:
    """
    mkdir -p baseATpositions_for baseATpositions_rev

    awk '{print "echo "\$5">baseATpositions_for/"\$1"_"\$2"_a"}' $maf_for | sh
    awk '{print "echo "\$1"\\t"\$2"\\t"\$7">baseATpositions_for/"\$1"_"\$2"_b"}' $maf_for | sh

    awk '{print "echo "\$5">baseATpositions_rev/"\$1"_"\$2"_a"}' $maf_rev | sh
    awk '{print "echo "\$1"\\t"\$2"\\t"\$7">baseATpositions_rev/"\$1"_"\$2"_b"}' $maf_rev | sh
    """
}

process mutate_and_merge {
    publishDir "data/base_information"

    input:
    path baseATpos_for
    path baseATpos_rev

    output:
    path "baseATpositions_for/merged_mutated_for", type: 'dir'
    path "baseATpositions_rev/merged_mutated_rev", type: 'dir'

    script:
    """
    mkdir -p baseATpositions_for/a_mutated baseATpositions_rev/a_mutated
    mkdir -p baseATpositions_for/b_mutated baseATpositions_rev/b_mutated
    mkdir -p baseATpositions_for/merged_mutated_for baseATpositions_rev/merged_mutated_rev

    for file in baseATpositions_for/*a; do
        name=\$(basename "\$file")
        awk -v FS="" '{for (i=1;i<=NF;i++) printf \$i"\\n"}' "\$file" |
        awk '{if (\$0 != "+" && \$0 != "-" && \$0 !~ /^[0-9]*\$/) print \$0}' \
        > baseATpositions_for/a_mutated/\${name}_mutated
    done

    for file in baseATpositions_rev/*a; do
        name=\$(basename "\$file")
        awk -v FS="" '{for (i=1;i<=NF;i++) printf \$i"\\n"}' "\$file" |
        awk '{if (\$0 != "+" && \$0 != "-" && \$0 !~ /^[0-9]*\$/) print \$0}' \
        > baseATpositions_rev/a_mutated/\${name}_mutated
    done

    for file in baseATpositions_for/*b; do
        name=\$(basename "\$file")
        awk '{print \$3}' "\$file" | sed 's/,/\\n/g' \
        > baseATpositions_for/b_mutated/\${name}_mutated
    done

    for file in baseATpositions_rev/*b; do
        name=\$(basename "\$file")
        awk '{print \$3}' "\$file" | sed 's/,/\\n/g' \
        > baseATpositions_rev/b_mutated/\${name}_mutated
    done

    for file in baseATpositions_for/*a; do
        id=\$(basename "\$file" | awk -F'_' '{print \$1"_"\$2"_"\$3}')
        paste baseATpositions_for/a_mutated/\${id}_a_mutated \
            baseATpositions_for/b_mutated/\${id}_b_mutated \
        > baseATpositions_for/merged_mutated_for/\${id}_merged
    done

    for file in baseATpositions_rev/*a; do
        id=\$(basename "\$file" | awk -F'_' '{print \$1"_"\$2"_"\$3}')
        paste baseATpositions_rev/a_mutated/\${id}_a_mutated \
            baseATpositions_rev/b_mutated/\${id}_b_mutated \
        > baseATpositions_rev/merged_mutated_rev/\${id}_merged
    done
    """
}

process umi_combining {
    publishDir "data/base_information"

    input:
    path bam
    path ribo_bed
    path merged_for
    path merged_rev
    path umi_combining_for_rscript
    path umi_combining_rev_rscript

    output:
    path "B2_deoxy_UMI_bases_all_for.tsv"
    path "B2_deoxy_UMI_bases_all_rev.tsv"
    path "B2_ribo_UMI_bases_all_for.tsv"
    path "B2_ribo_UMI_bases_all_rev.tsv"

    script:
    """
    samtools view $bam |cut -f1,17,19 > channel_read_number_all 

    for r in \$(cut -f1 $ribo_bed | sort -u); do
      Rscript $umi_combining_for_rscript \$r ./
      Rscript $umi_combining_rev_rscript \$r ./
    done
    """
}

process umi_pairing {
    publishDir "data/"

    input:
    path umi_for_ctl
    path umi_rev_ctl
    path umi_for_ribo
    path umi_rev_ribo
    path umi_pairing_rscript

    output:
    path "ctl_paired_output.tsv"
    path "ribo_paired_output.tsv"

    script:
    """
    Rscript $umi_pairing_rscript $umi_for_ctl $umi_rev_ctl ctl_paired_output.tsv
    Rscript $umi_pairing_rscript $umi_for_ribo $umi_rev_ribo ribo_paired_output.tsv
    """
}

process summary_stats {
    publishDir "data/"

    input:
    path ctl_pairs
    path ribo_pairs
    path ctl_all_for
    path ctl_all_rev
    path ribo_all_for
    path ribo_all_rev

    output:
    path "summary_ctl.txt"
    path "summary_ribo.txt"

    script:
    """
    echo "ctl pairs: \$(grep -v "*" $ctl_pairs | grep -v rname | wc -l)" > summary_ctl.txt
    echo "ctl forward total: \$(grep -v "*" $ctl_all_for | grep -v rname | wc -l)" >> summary_ctl.txt
    echo "ctl reverse total: \$(grep -v "*" $ctl_all_rev | grep -v rname | wc -l)" >> summary_ctl.txt

    echo "ribo pairs: \$(grep -v "*" $ribo_pairs | grep -v rname | wc -l)" > summary_ribo.txt
    echo "ribo forward total: \$(grep -v "*" $ribo_all_for | grep -v rname | wc -l)" >> summary_ribo.txt
    echo "ribo reverse total: \$(grep -v "*" $ribo_all_rev | grep -v rname | wc -l)" >> summary_ribo.txt
    """
}

process pileup_trinu {
    publishDir "data/"

    input:
    path bam_for
    path bam_rev
    path ref
    path ribo_bed

    output:
    path "mpileup_output_for_ribo.maf"
    path "mpileup_output_rev_ribo.maf"

    script:
    """
    samtools mpileup $bam_for -f $ref --no-output-ins --no-output-del --no-output-ends \
      -l $ribo_bed -o mpileup_output_for_ribo.maf --output-QNAME -A -B --max-depth 0

    samtools mpileup $bam_rev -f $ref --no-output-ins --no-output-del --no-output-ends \
      -l $ribo_bed -o mpileup_output_rev_ribo.maf --output-QNAME -A -B --max-depth 0
    """
}

process split_bases_ribo {
    publishDir "data/base_information_ribo"

    input:
    path maf_for
    path maf_rev

    output:
    path "baseATpositions_for", type: 'dir'
    path "baseATpositions_rev", type: 'dir'

    script:
    """
    mkdir -p baseATpositions_for baseATpositions_rev

    awk '{print "echo "\$5">baseATpositions_for/"\$1"_"\$2"_a"}' $maf_for | sh
    awk '{print "echo "\$1"\\t"\$2"\\t"\$7">baseATpositions_for/"\$1"_"\$2"_b"}' $maf_for | sh

    awk '{print "echo "\$5">baseATpositions_rev/"\$1"_"\$2"_a"}' $maf_rev | sh
    awk '{print "echo "\$1"\\t"\$2"\\t"\$7">baseATpositions_rev/"\$1"_"\$2"_b"}' $maf_rev | sh
    """
}

process mutate_and_merge_ribo {
    publishDir "data/base_information_ribo"

    input:
    path baseATpos_for
    path baseATpos_rev

    output:
    path "baseATpositions_for/merged_mutated_for", type: 'dir'
    path "baseATpositions_rev/merged_mutated_rev", type: 'dir'

    script:
    """
    mkdir -p baseATpositions_for/a_mutated baseATpositions_rev/a_mutated
    mkdir -p baseATpositions_for/b_mutated baseATpositions_rev/b_mutated
    mkdir -p baseATpositions_for/merged_mutated_for baseATpositions_rev/merged_mutated_rev

    for file in baseATpositions_for/*a; do
        name=\$(basename "\$file")
        awk -v FS="" '{for (i=1;i<=NF;i++) printf \$i"\\n"}' "\$file" |
        awk '{if (\$0 != "+" && \$0 != "-" && \$0 !~ /^[0-9]*\$/) print \$0}' \
        > baseATpositions_for/a_mutated/\${name}_mutated
    done

    for file in baseATpositions_rev/*a; do
        name=\$(basename "\$file")
        awk -v FS="" '{for (i=1;i<=NF;i++) printf \$i"\\n"}' "\$file" |
        awk '{if (\$0 != "+" && \$0 != "-" && \$0 !~ /^[0-9]*\$/) print \$0}' \
        > baseATpositions_rev/a_mutated/\${name}_mutated
    done

    for file in baseATpositions_for/*b; do
        name=\$(basename "\$file")
        awk '{print \$3}' "\$file" | sed 's/,/\\n/g' \
        > baseATpositions_for/b_mutated/\${name}_mutated
    done

    for file in baseATpositions_rev/*b; do
        name=\$(basename "\$file")
        awk '{print \$3}' "\$file" | sed 's/,/\\n/g' \
        > baseATpositions_rev/b_mutated/\${name}_mutated
    done

    for file in baseATpositions_for/*a; do
        id=\$(basename "\$file" | awk -F'_' '{print \$1"_"\$2"_"\$3}')
        paste baseATpositions_for/a_mutated/\${id}_a_mutated \
            baseATpositions_for/b_mutated/\${id}_b_mutated \
        > baseATpositions_for/merged_mutated_for/\${id}_merged
    done

    for file in baseATpositions_rev/*a; do
        id=\$(basename "\$file" | awk -F'_' '{print \$1"_"\$2"_"\$3}')
        paste baseATpositions_rev/a_mutated/\${id}_a_mutated \
            baseATpositions_rev/b_mutated/\${id}_b_mutated \
        > baseATpositions_rev/merged_mutated_rev/\${id}_merged
    done
    """
}

process calculate_trinu_acc {
    publishDir "data/trinu_acc"

    input:
    path trinu_base_for
    path trinu_base_rev

    output:
    path "merged_mutated_*/B2_*_merged_acc.tsv"

    script:
    """
    for l in $trinu_base_for/B2*; do
        out="\${l}_acc.tsv"
        awk '{if (\$1 == ".") sum++; count++} END {print sum"\\t"count"\\t"sum/count"\\t"count - sum}' "\$l" > "\$out"
    done

    for l in $trinu_base_rev/B2*; do
        out="\${l}_acc.tsv"
        awk '{if (\$1 == ",") sum++; count++} END {print sum"\\t"count"\\t"sum/count"\\t"count - sum}' "\$l" > "\$out"
    done
    """
}

workflow {
    // Define channels for inputs
    bam_ch       = Channel.fromPath(params.bam)
    ref_ch       = Channel.fromPath(params.ref)
    umi_ch       = Channel.fromPath(params.umi_bed)
    ribo_ch      = Channel.fromPath(params.ribo_bed)
    umi_r_for_ch = Channel.fromPath(params.umi_r_for)
    umi_r_rev_ch = Channel.fromPath(params.umi_r_rev)
    r_pairing_ch = Channel.fromPath(params.umi_r_pairing)    

    // Stepwise pipeline
    fastq_ch     = bam2fastq(bam_ch)
    aln_ch       = minimap2_align(fastq_ch, ref_ch)
    primary_for_ch   = filter_primary_for(aln_ch)
    primary_rev_ch   = filter_primary_rev(aln_ch)
    conca_ch     = filter_concatemers(primary_for_ch, primary_rev_ch)
    essential_ch = essential_regions(primary_for_ch, primary_rev_ch)
    inter_ch     = intersect_reads(conca_ch, essential_ch)
    bam_filtered_ch      = filter_bam(primary_for_ch, primary_rev_ch, inter_ch)
    index_ch     = index_filtered_bam(bam_filtered_ch)
    maf_ch     = pileup(bam_filtered_ch, ref_ch, umi_ch)
    bases_ch   = split_bases(maf_ch)
    merged_ch  = mutate_and_merge(bases_ch)
    umi_ch     = umi_combining(bam_ch, ribo_ch, merged_ch, umi_r_for_ch, umi_r_rev_ch)
    pair_ch    = umi_pairing(umi_ch, r_pairing_ch)
    summary_stats(pair_ch, umi_ch)
    maf_ribo_ch   = pileup_trinu(bam_filtered_ch, ref_ch, ribo_ch)
    bases_ribo_ch = split_bases_ribo(maf_ribo_ch)
    merged_ribo_ch  = mutate_and_merge_ribo(bases_ribo_ch)
    calculate_trinu_acc(merged_ribo_ch)
}
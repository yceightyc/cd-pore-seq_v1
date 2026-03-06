#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.bam_r10   = "./input/R10_simplex.bam"
params.bam_hd    = "./input/HD_simplex.bam"
params.ref       = "./input/reference.fasta"
params.roi_bed   = "./input/barcode_ROI_UMI.bed"
params.fil_barcode_rscript = "./filter_barcode_calls.R"
params.umi_r_for = "./UMI_permut_combining_bases_for.R"
params.umi_r_rev = "./UMI_permut_combining_bases_rev.R"
params.bandu_r_for = "./pull_barcoded_filtered_for.R"
params.bandu_r_rev = "./pull_barcoded_filtered_rev.R"
params.pairing_r = "./pairing_up.R"
params.degen_r   = "./UMI_degeneracy_analysis.R"
params.plot_r    = "./plot_degeneracy.R"

process bam2fastq {
    publishDir "data/"

    input:
    path bam_r10
    path bam_hd

    output:
    path "R10_simplex.fastq"
    path "HD_simplex.fastq"

    script:
    """
    samtools fastq -T "*" $bam_r10 > R10_simplex.fastq
    samtools fastq -T "*" $bam_hd > HD_simplex.fastq
    """
}

process minimap2_align {
    publishDir "data/"

    input:
    path fastq_r10
    path fastq_hd
    path ref

    output:
    path "minimap2_k5w3_R10.bam"
    path "minimap2_k5w3_HD.bam"

    script:
    """
    minimap2 --MD -y -a $ref $fastq_r10 -k 5 -w 3 | \
      samtools sort -T tmp -o minimap2_k5w3_R10.bam
    
    minimap2 --MD -y -a $ref $fastq_hd -k 5 -w 3 | \
      samtools sort -T tmp -o minimap2_k5w3_HD.bam
    """
}

process filter_primary_for {
    publishDir "data/"

    input:
    path bam_r10
    path bam_hd

    output:
    path "minimap2_pri_forward_R10.bam"
    path "minimap2_pri_forward_HD.bam"

    script:
    """
    samtools view -bS -F 256 -F 16 -F 4 -F 2048 $bam_r10 | \
      samtools sort -T tMp -o minimap2_pri_forward_R10.bam ;
    
    samtools view -bS -F 256 -F 16 -F 4 -F 2048 $bam_hd | \
      samtools sort -T tMp -o minimap2_pri_forward_HD.bam ;
    """
}

process filter_primary_rev {
    publishDir "data/"

    input:
    path bam_r10
    path bam_hd

    output:
    path "minimap2_pri_reverse_R10.bam"
    path "minimap2_pri_reverse_HD.bam"

    script:
    """
    samtools view -bS -F 256 -f 16 -F 4 -F 2048 $bam_r10 | \
      samtools sort -T Tmp -o minimap2_pri_reverse_R10.bam ;

    samtools view -bS -F 256 -f 16 -F 4 -F 2048 $bam_hd | \
      samtools sort -T Tmp -o minimap2_pri_reverse_HD.bam ;
    """
}

process filter_concatemers {
    publishDir "data/"

    input:
    path bam_r10_for
    path bam_hd_for
    path bam_r10_rev
    path bam_hd_rev

    output:
    path "noconca_RN_for_R10.rn"
    path "noconca_RN_for_HD.rn"
    path "noconca_RN_rev_R10.rn"
    path "noconca_RN_rev_HD.rn"

 
    script:
    """
    # The correct full size of a read including adapters should be 123 + 36 + 27 = 186, set threshold as 190 to allow small insertions
    samtools view $bam_r10_for | \
      awk '{print \$1"\\t"\$4"\\t"length(\$10)}' | \
      awk '{print \$1"\\t"\$2+\$3}' | \
      awk '{if (\$2<=190) print \$0}' | \
      cut -f1 > noconca_RN_for_R10.rn

    samtools view $bam_hd_for | \
      awk '{print \$1"\\t"\$4"\\t"length(\$10)}' | \
      awk '{print \$1"\\t"\$2+\$3}' | \
      awk '{if (\$2<=190) print \$0}' | \
      cut -f1 > noconca_RN_for_HD.rn
    
    samtools view $bam_r10_rev | \
      awk '{print \$1"\\t"\$4"\\t"length(\$10)}' | \
      awk '{print \$1"\\t"\$2+\$3}' | \
      awk '{if (\$2<=190) print \$0}' | \
      cut -f1 > noconca_RN_rev_R10.rn
    
    samtools view $bam_hd_rev | \
      awk '{print \$1"\\t"\$4"\\t"length(\$10)}' | \
      awk '{print \$1"\\t"\$2+\$3}' | \
      awk '{if (\$2<=190) print \$0}' | \
      cut -f1 > noconca_RN_rev_HD.rn
    """
}

process essential_regions {
    publishDir "data/"

    input:
    path bam_r10_for
    path bam_hd_for
    path bam_r10_rev
    path bam_hd_rev

    output:
    path "essential_RN_for_R10.rn"
    path "essential_RN_for_HD.rn"
    path "essential_RN_rev_R10.rn"
    path "essential_RN_rev_HD.rn"

    script:
    """
    bedtools bamtobed -i $bam_r10_for | \
      awk '{if (\$2 <= 8 && \$3 >= 94) print \$0}' | \
      cut -f4 > essential_RN_for_R10.rn

    bedtools bamtobed -i $bam_hd_for | \
      awk '{if (\$2 <= 8 && \$3 >= 94) print \$0}' | \
      cut -f4 > essential_RN_for_HD.rn
    
    bedtools bamtobed -i $bam_r10_rev | \
      awk '{if (\$2 <= 8 && \$3 >= 94) print \$0}' | \
      cut -f4 > essential_RN_rev_R10.rn

    bedtools bamtobed -i $bam_hd_rev | \
      awk '{if (\$2 <= 8 && \$3 >= 94) print \$0}' | \
      cut -f4 > essential_RN_rev_HD.rn
    """
}

process intersect_reads {
    publishDir "data/"

    input:
    path rn1_for_r10
    path rn1_for_hd
    path rn1_rev_r10
    path rn1_rev_hd
    path rn2_for_r10
    path rn2_for_hd
    path rn2_rev_r10
    path rn2_rev_hd

    output:
    path "non_conca_with_essential_RN_for_R10.rn"
    path "non_conca_with_essential_RN_for_HD.rn"
    path "non_conca_with_essential_RN_rev_R10.rn"
    path "non_conca_with_essential_RN_rev_HD.rn"

    script:
    """
    cat $rn1_for_r10 $rn2_for_r10 | sort | uniq -d > non_conca_with_essential_RN_for_R10.rn
    cat $rn1_for_hd $rn2_for_hd | sort | uniq -d > non_conca_with_essential_RN_for_HD.rn
    cat $rn1_rev_r10 $rn2_rev_r10 | sort | uniq -d > non_conca_with_essential_RN_rev_R10.rn
    cat $rn1_rev_hd $rn2_rev_hd | sort | uniq -d > non_conca_with_essential_RN_rev_HD.rn
    """
}

process self_annealing {
    publishDir "data/"

    input:
    path bam_r10_for
    path bam_hd_for
    path bam_r10_rev
    path bam_hd_rev

    output:
    path "selfannealing_RN_for_R10.rn"
    path "selfannealing_RN_for_HD.rn"
    path "selfannealing_RN_rev_R10.rn"
    path "selfannealing_RN_rev_HD.rn"

    script:
    """
    samtools view $bam_r10_for | \
      grep "GCGTGACATCGATGACACGC" | \
      cut -f1 > selfannealing_RN_for_R10.rn

    samtools view $bam_hd_for | \
      grep "GCGTGACATCGATGACACGC" | \
      cut -f1 > selfannealing_RN_for_HD.rn
    
    samtools view $bam_r10_rev | \
      grep "GCGTGACATCGATGACACGC" | \
      cut -f1 > selfannealing_RN_rev_R10.rn

    samtools view $bam_hd_rev | \
      grep "GCGTGACATCGATGACACGC" | \
      cut -f1 > selfannealing_RN_rev_HD.rn
    """
}

process remove_sa {
    publishDir "data/"

    input:
    path sa_rn_r10_for
    path sa_rn_hd_for
    path sa_rn_r10_rev
    path sa_rn_hd_rev
    path rn_for_r10
    path rn_for_hd
    path rn_rev_r10
    path rn_rev_hd

    output:
    path "non_conca_with_essential_no_selfannealing_RN_for_R10.rn"
    path "non_conca_with_essential_no_selfannealing_RN_for_HD.rn"
    path "non_conca_with_essential_no_selfannealing_RN_rev_R10.rn"
    path "non_conca_with_essential_no_selfannealing_RN_rev_HD.rn"

    script:
    """
    cat $sa_rn_r10_for $sa_rn_r10_for $rn_for_r10|sort| \
      uniq -u > non_conca_with_essential_no_selfannealing_RN_for_R10.rn
    cat $sa_rn_hd_for $sa_rn_hd_for $rn_for_hd|sort| \
      uniq -u > non_conca_with_essential_no_selfannealing_RN_for_HD.rn
    cat $sa_rn_r10_rev $sa_rn_r10_rev $rn_rev_r10|sort| \
      uniq -u > non_conca_with_essential_no_selfannealing_RN_rev_R10.rn
    cat $sa_rn_hd_rev $sa_rn_hd_rev $rn_rev_hd|sort| \
      uniq -u > non_conca_with_essential_no_selfannealing_RN_rev_HD.rn
    """
}

process filter_bam {
    publishDir "data/"

    input:
    path bam_for_r10
    path bam_for_hd
    path bam_rev_r10
    path bam_rev_hd
    path rn_for_r10
    path rn_for_hd
    path rn_rev_r10
    path rn_rev_hd

    output:
    path "filterd_for_R10.bam"
    path "filterd_for_HD.bam"
    path "filterd_rev_R10.bam"
    path "filterd_rev_HD.bam"
    
    script:
    """
    samtools view -N $rn_for_r10 -o filterd_for_R10.bam $bam_for_r10
    samtools view -N $rn_for_hd -o filterd_for_HD.bam $bam_for_hd
    samtools view -N $rn_rev_r10 -o filterd_rev_R10.bam $bam_rev_r10
    samtools view -N $rn_rev_hd -o filterd_rev_HD.bam $bam_rev_hd
    """
}

process index_filtered_bam {
    publishDir "data/"

    input:
    path bam_for_r10
    path bam_for_hd
    path bam_rev_r10
    path bam_rev_hd

    output:
    path "filterd_for_R10.bam.bai"
    path "filterd_for_HD.bam.bai"
    path "filterd_rev_R10.bam.bai"
    path "filterd_rev_HD.bam.bai"

    script:
    """
    samtools index $bam_for_r10
    samtools index $bam_for_hd
    samtools index $bam_rev_r10
    samtools index $bam_rev_hd
    """
}

process pileup {
    publishDir "data/"

    input:
    path bam_for_r10
    path bam_for_hd
    path bam_rev_r10
    path bam_rev_hd
    path ref
    path roi_bed

    output:
    path "mpileup_output_for_R10.maf"
    path "mpileup_output_for_HD.maf"
    path "mpileup_output_rev_R10.maf"
    path "mpileup_output_rev_HD.maf"

    script:
    """
    samtools mpileup $bam_for_r10 -f $ref --no-output-ins --no-output-del --no-output-ends \
      -l $roi_bed -o mpileup_output_for_R10.maf --output-QNAME -A -B --max-depth 0

    samtools mpileup $bam_for_hd -f $ref --no-output-ins --no-output-del --no-output-ends \
      -l $roi_bed -o mpileup_output_for_HD.maf --output-QNAME -A -B --max-depth 0

    samtools mpileup $bam_rev_r10 -f $ref --no-output-ins --no-output-del --no-output-ends \
      -l $roi_bed -o mpileup_output_rev_R10.maf --output-QNAME -A -B --max-depth 0

    samtools mpileup $bam_rev_hd -f $ref --no-output-ins --no-output-del --no-output-ends \
      -l $roi_bed -o mpileup_output_rev_HD.maf --output-QNAME -A -B --max-depth 0
    """
}

process split_bases {
    publishDir "data/base_information"

    input:
    path maf_for_r10
    path maf_for_hd
    path maf_rev_r10
    path maf_rev_hd

    output:
    path "baseATpositions_for_R10", type: 'dir'
    path "baseATpositions_for_HD", type: 'dir'
    path "baseATpositions_rev_R10", type: 'dir'
    path "baseATpositions_rev_HD", type: 'dir'

    script:
    """
    mkdir -p baseATpositions_for_R10 baseATpositions_for_HD
    mkdir -p baseATpositions_rev_R10 baseATpositions_rev_HD

    awk '{print "echo "\$5">baseATpositions_for_R10/"\$1"_"\$2"_a"}' $maf_for_r10 | sh
    awk '{print "echo "\$1"\\t"\$2"\\t"\$7">baseATpositions_for_R10/"\$1"_"\$2"_b"}' $maf_for_r10 | sh

    awk '{print "echo "\$5">baseATpositions_for_HD/"\$1"_"\$2"_a"}' $maf_for_hd | sh
    awk '{print "echo "\$1"\\t"\$2"\\t"\$7">baseATpositions_for_HD/"\$1"_"\$2"_b"}' $maf_for_hd | sh

    awk '{print "echo "\$5">baseATpositions_rev_R10/"\$1"_"\$2"_a"}' $maf_rev_r10 | sh
    awk '{print "echo "\$1"\\t"\$2"\\t"\$7">baseATpositions_rev_R10/"\$1"_"\$2"_b"}' $maf_rev_r10 | sh

    awk '{print "echo "\$5">baseATpositions_rev_HD/"\$1"_"\$2"_a"}' $maf_rev_hd | sh
    awk '{print "echo "\$1"\\t"\$2"\\t"\$7">baseATpositions_rev_HD/"\$1"_"\$2"_b"}' $maf_rev_hd | sh
    """
}

process mutate_and_merge {
    publishDir "data/base_information"

    input:
    path baseATpos_for_r10
    path baseATpos_for_hd
    path baseATpos_rev_r10
    path baseATpos_rev_hd

    output:
    path "baseATpositions_for_R10/merged_mutated_for_R10", type: 'dir'
    path "baseATpositions_for_HD/merged_mutated_for_HD", type: 'dir'
    path "baseATpositions_rev_R10/merged_mutated_rev_R10", type: 'dir'
    path "baseATpositions_rev_HD/merged_mutated_rev_HD", type: 'dir'

    script:
    """
    mkdir -p baseATpositions_for_R10/a_mutated baseATpositions_rev_R10/a_mutated
    mkdir -p baseATpositions_for_R10/b_mutated baseATpositions_rev_R10/b_mutated
    mkdir -p baseATpositions_for_R10/merged_mutated_for_R10 baseATpositions_rev_R10/merged_mutated_rev_R10

    mkdir -p baseATpositions_for_HD/a_mutated baseATpositions_rev_HD/a_mutated
    mkdir -p baseATpositions_for_HD/b_mutated baseATpositions_rev_HD/b_mutated
    mkdir -p baseATpositions_for_HD/merged_mutated_for_HD baseATpositions_rev_HD/merged_mutated_rev_HD

    for file in baseATpositions_for_R10/*a; do
        name=\$(basename "\$file")
        awk -v FS="" '{for (i=1;i<=NF;i++) printf \$i"\\n"}' "\$file" |
        awk '{if (\$0 != "+" && \$0 != "-" && \$0 !~ /^[0-9]*\$/) print \$0}' \
        > baseATpositions_for_R10/a_mutated/\${name}_mutated
    done

    for file in baseATpositions_rev_R10/*a; do
        name=\$(basename "\$file")
        awk -v FS="" '{for (i=1;i<=NF;i++) printf \$i"\\n"}' "\$file" |
        awk '{if (\$0 != "+" && \$0 != "-" && \$0 !~ /^[0-9]*\$/) print \$0}' \
        > baseATpositions_rev_R10/a_mutated/\${name}_mutated
    done

    for file in baseATpositions_for_R10/*b; do
        name=\$(basename "\$file")
        awk '{print \$3}' "\$file" | sed 's/,/\\n/g' \
        > baseATpositions_for_R10/b_mutated/\${name}_mutated
    done

    for file in baseATpositions_rev_R10/*b; do
        name=\$(basename "\$file")
        awk '{print \$3}' "\$file" | sed 's/,/\\n/g' \
        > baseATpositions_rev_R10/b_mutated/\${name}_mutated
    done

    for file in baseATpositions_for_R10/*a; do
        id=\$(basename "\$file" | awk -F'_' '{print \$1"_"\$2}')
        paste baseATpositions_for_R10/a_mutated/\${id}_a_mutated \
            baseATpositions_for_R10/b_mutated/\${id}_b_mutated \
        > baseATpositions_for_R10/merged_mutated_for_R10/\${id}_merged
    done

    for file in baseATpositions_rev_R10/*a; do
        id=\$(basename "\$file" | awk -F'_' '{print \$1"_"\$2}')
        paste baseATpositions_rev_R10/a_mutated/\${id}_a_mutated \
            baseATpositions_rev_R10/b_mutated/\${id}_b_mutated \
        > baseATpositions_rev_R10/merged_mutated_rev_R10/\${id}_merged
    done

        for file in baseATpositions_for_HD/*a; do
        name=\$(basename "\$file")
        awk -v FS="" '{for (i=1;i<=NF;i++) printf \$i"\\n"}' "\$file" |
        awk '{if (\$0 != "+" && \$0 != "-" && \$0 !~ /^[0-9]*\$/) print \$0}' \
        > baseATpositions_for_HD/a_mutated/\${name}_mutated
    done

    for file in baseATpositions_rev_HD/*a; do
        name=\$(basename "\$file")
        awk -v FS="" '{for (i=1;i<=NF;i++) printf \$i"\\n"}' "\$file" |
        awk '{if (\$0 != "+" && \$0 != "-" && \$0 !~ /^[0-9]*\$/) print \$0}' \
        > baseATpositions_rev_HD/a_mutated/\${name}_mutated
    done

    for file in baseATpositions_for_HD/*b; do
        name=\$(basename "\$file")
        awk '{print \$3}' "\$file" | sed 's/,/\\n/g' \
        > baseATpositions_for_HD/b_mutated/\${name}_mutated
    done

    for file in baseATpositions_rev_HD/*b; do
        name=\$(basename "\$file")
        awk '{print \$3}' "\$file" | sed 's/,/\\n/g' \
        > baseATpositions_rev_HD/b_mutated/\${name}_mutated
    done

    for file in baseATpositions_for_HD/*a; do
        id=\$(basename "\$file" | awk -F'_' '{print \$1"_"\$2}')
        paste baseATpositions_for_HD/a_mutated/\${id}_a_mutated \
            baseATpositions_for_HD/b_mutated/\${id}_b_mutated \
        > baseATpositions_for_HD/merged_mutated_for_HD/\${id}_merged
    done

    for file in baseATpositions_rev_HD/*a; do
        id=\$(basename "\$file" | awk -F'_' '{print \$1"_"\$2}')
        paste baseATpositions_rev_HD/a_mutated/\${id}_a_mutated \
            baseATpositions_rev_HD/b_mutated/\${id}_b_mutated \
        > baseATpositions_rev_HD/merged_mutated_rev_HD/\${id}_merged
    done
    """
}

process barcode_integrity {
    publishDir "data/"

    input:
    path dir_merged_mutated_for_r10
    path dir_merged_mutated_for_hd
    path dir_merged_mutated_rev_r10
    path dir_merged_mutated_rev_hd
    path fil_barcode_rscript
    path roi_bed

    output:
    path "barcode_filtered_for_R10.tsv"
    path "barcode_filtered_for_HD.tsv"
    path "barcode_filtered_rev_R10.tsv"
    path "barcode_filtered_rev_HD.tsv"

    script:
    """
    less $roi_bed |cut -f1|sort -u|\
      awk '{print "Rscript $fil_barcode_rscript $dir_merged_mutated_for_r10 "\$0}'|sh ;
    cat $dir_merged_mutated_for_r10/*_filtered > barcode_filtered_for_R10.tsv

    less $roi_bed |cut -f1|sort -u|\
      awk '{print "Rscript $fil_barcode_rscript $dir_merged_mutated_for_hd "\$0}'|sh ;
    cat $dir_merged_mutated_for_hd/*_filtered > barcode_filtered_for_HD.tsv

    less $roi_bed |cut -f1|sort -u|\
      awk '{print "Rscript $fil_barcode_rscript $dir_merged_mutated_rev_r10 "\$0}'|sh ;
    cat $dir_merged_mutated_rev_r10/*_filtered > barcode_filtered_rev_R10.tsv
  
    less $roi_bed |cut -f1|sort -u|\
      awk '{print "Rscript $fil_barcode_rscript $dir_merged_mutated_rev_hd "\$0}'|sh ;
    cat $dir_merged_mutated_rev_hd/*_filtered > barcode_filtered_rev_HD.tsv
    """
}

process umi_combining {
    publishDir "data/"

    input:
    path bam_r10
    path bam_hd
    path roi_bed
    path dir_merged_mutated_for_r10
    path dir_merged_mutated_for_hd
    path dir_merged_mutated_rev_r10
    path dir_merged_mutated_rev_hd
    path umi_combining_for_rscript
    path umi_combining_rev_rscript

    output:
    path "UMI_bases_aggregated_for_R10.tsv"
    path "UMI_bases_aggregated_for_HD.tsv"
    path "UMI_bases_aggregated_rev_R10.tsv"
    path "UMI_bases_aggregated_rev_HD.tsv"

    script:
    """
    samtools view $bam_r10 |cut -f1,17,19 > $dir_merged_mutated_for_r10/channel_read_number_all 
    samtools view $bam_hd |cut -f1,17,19 > $dir_merged_mutated_for_hd/channel_read_number_all
    cp $dir_merged_mutated_for_r10/channel_read_number_all $dir_merged_mutated_rev_r10/channel_read_number_all
    cp $dir_merged_mutated_for_hd/channel_read_number_all $dir_merged_mutated_rev_hd/channel_read_number_all
    
    # Count *_all.tsv files in each folder

    if [ "\$(ls $dir_merged_mutated_for_r10/*_all.tsv|wc -l)" -eq 64 ] && [ "\$(ls $dir_merged_mutated_for_hd/*_all.tsv|wc -l)" -eq 64 ] && \
       [ "\$(ls $dir_merged_mutated_rev_r10/*_all.tsv|wc -l)" -eq 64 ] && [ "\$(ls $dir_merged_mutated_rev_hd/*_all.tsv|wc -l)" -eq 64 ]; then
        : 
    else
        for r in \$(cut -f1 $roi_bed | sort -u); do
          Rscript $umi_combining_for_rscript \$r $dir_merged_mutated_for_r10
          Rscript $umi_combining_for_rscript \$r $dir_merged_mutated_for_hd
          Rscript $umi_combining_rev_rscript \$r $dir_merged_mutated_rev_r10
          Rscript $umi_combining_rev_rscript \$r $dir_merged_mutated_rev_hd
        done
    fi

    ls $dir_merged_mutated_for_r10/*_all.tsv | \
      awk -F'[/|_]' '{print \$0"\t"\$5}'| \
      while read a b;do awk '{print "'\$b'\t"\$0}' \$a;done| \
      grep -wv readname > UMI_bases_aggregated_for_R10.tsv

    ls $dir_merged_mutated_for_hd/*_all.tsv | \
      awk -F'[/|_]' '{print \$0"\t"\$5}'| \
      while read a b;do awk '{print "'\$b'\t"\$0}' \$a;done| \
      grep -wv readname > UMI_bases_aggregated_for_HD.tsv
    
    ls $dir_merged_mutated_rev_r10/*_all.tsv | \
      awk -F'[/|_]' '{print \$0"\t"\$5}'| \
      while read a b;do awk '{print "'\$b'\t"\$0}' \$a;done| \
      grep -wv readname > UMI_bases_aggregated_rev_R10.tsv

    ls $dir_merged_mutated_rev_hd/*_all.tsv | \
      awk -F'[/|_]' '{print \$0"\t"\$5}'| \
      while read a b;do awk '{print "'\$b'\t"\$0}' \$a;done| \
      grep -wv readname > UMI_bases_aggregated_rev_HD.tsv
    """
}

process final_filter {
    publishDir "data/"

    input:
    path barcode_filtered_for_rn_r10
    path barcode_filtered_for_rn_hd
    path barcode_filtered_rev_rn_r10
    path barcode_filtered_rev_rn_hd
    path umi_base_for_r10
    path umi_base_for_hd
    path umi_base_rev_r10
    path umi_base_rev_hd
    path bandu_for_rscript
    path bandu_rev_rscript
    
    output:
    path "final_UMI_for_R10.tsv"
    path "final_UMI_for_HD.tsv"
    path "final_UMI_rev_R10.tsv"
    path "final_UMI_rev_HD.tsv"

    script:
    """
    Rscript $bandu_for_rscript $barcode_filtered_for_rn_r10 $umi_base_for_r10 final_UMI_for_R10.tsv
    Rscript $bandu_for_rscript $barcode_filtered_for_rn_hd $umi_base_for_hd final_UMI_for_HD.tsv
    Rscript $bandu_rev_rscript $barcode_filtered_rev_rn_r10 $umi_base_rev_r10 final_UMI_rev_R10.tsv
    Rscript $bandu_rev_rscript $barcode_filtered_rev_rn_hd $umi_base_rev_hd final_UMI_rev_HD.tsv
    """
}

process umi_pairing {
    publishDir "data/"

    input:
    path umi_for_r10
    path umi_for_hd
    path umi_rev_r10
    path umi_rev_hd
    path umi_pairing_rscript

    output:
    path "paired_output_R10.tsv"
    path "paired_output_HD.tsv"

    script:
    """
    Rscript $umi_pairing_rscript $umi_for_r10 $umi_rev_r10 paired_output_R10.tsv
    Rscript $umi_pairing_rscript $umi_for_hd $umi_rev_hd paired_output_HD.tsv
    """
}

process umi_degeneracy {
    publishDir "data/" , mode:"copy"

    input:
    path r10_all_for
    path hd_all_for
    path r10_all_rev
    path hd_all_rev
    path umi_deg_rscript

    output:
    path "UMIdegeneracy_R10.tsv"
    path "UMIdegeneracy_HD.tsv"

    script:
    """
    Rscript $umi_deg_rscript $r10_all_for $r10_all_rev UMIdegeneracy_R10.tsv
    Rscript $umi_deg_rscript $hd_all_for $hd_all_rev UMIdegeneracy_HD.tsv
    """
}

process plot_degeneracy {
    publishDir "data/" , mode:"copy"

    input:
    path degen_r10
    path degen_hd
    path plot_deg_rscript

    output:
    path "UMI_degeneracy_R10.pdf"
    path "UMI_degeneracy_HD.pdf"

    script:
    """
    Rscript $plot_deg_rscript $degen_r10 UMI_degeneracy_R10.pdf
    Rscript $plot_deg_rscript $degen_hd UMI_degeneracy_HD.pdf
    """
}

process summary_stats {
    publishDir "data/"

    input:
    path r10_pairs
    path hd_pairs
    path r10_all_for
    path hd_all_for
    path r10_all_rev
    path hd_all_rev

    output:
    path "summary_R10.txt"
    path "summary_HD.txt"

    script:
    """
    echo "R10 pairs: \$(grep -v "*" $r10_pairs |grep -v rname| wc -l)" > summary_R10.txt
    echo "R10 forward total: \$(grep -v "*" $r10_all_for | grep -v rname | wc -l)" >> summary_R10.txt
    echo "R10 reverse total: \$(grep -v "*" $r10_all_rev | grep -v rname | wc -l)" >> summary_R10.txt

    echo "HD pairs: \$(grep -v "*" $hd_pairs |grep -v rname| wc -l)" > summary_HD.txt
    echo "HD forward total: \$(grep -v "*" $hd_all_for | grep -v rname | wc -l)" >> summary_HD.txt
    echo "HD reverse total: \$(grep -v "*" $hd_all_rev | grep -v rname | wc -l)" >> summary_HD.txt
    """
}

workflow {
    // Define channels for inputs
    bam_r10_ch   = Channel.fromPath(params.bam_r10)
    bam_hd_ch    = Channel.fromPath(params.bam_hd)
    ref_ch       = Channel.fromPath(params.ref)   
    bed_ch       = Channel.fromPath(params.roi_bed)
    fil_bar_r_ch = Channel.fromPath(params.fil_barcode_rscript)
    umi_r_for_ch = Channel.fromPath(params.umi_r_for)
    umi_r_rev_ch = Channel.fromPath(params.umi_r_rev)
    bandu_r_for_ch = Channel.fromPath(params.bandu_r_for)
    bandu_r_rev_ch = Channel.fromPath(params.bandu_r_rev)
    pairing_r_ch = Channel.fromPath(params.pairing_r)
    degen_r_ch   = Channel.fromPath(params.degen_r)
    plot_r_ch    = Channel.fromPath(params.plot_r)

    // Stepwise pipeline
    fastq_ch     = bam2fastq(bam_r10_ch, bam_hd_ch)
    aln_ch       = minimap2_align(fastq_ch, ref_ch)
    primary_for_ch   = filter_primary_for(aln_ch)
    primary_rev_ch   = filter_primary_rev(aln_ch)
    conca_ch     = filter_concatemers(primary_for_ch, primary_rev_ch)
    essential_ch = essential_regions(primary_for_ch, primary_rev_ch)
    inter_ch     = intersect_reads(conca_ch, essential_ch)
    selfann_ch   = self_annealing(primary_for_ch, primary_rev_ch)
    rn_ch        = remove_sa(selfann_ch, inter_ch)
    bam_filtered_ch      = filter_bam(primary_for_ch, primary_rev_ch, rn_ch)
    index_ch     = index_filtered_bam(bam_filtered_ch)
    maf_ch       = pileup(bam_filtered_ch, ref_ch, bed_ch)
    bases_ch     = split_bases(maf_ch) 
    merged_ch    = mutate_and_merge(bases_ch)
    bar_fil_ch   = barcode_integrity(merged_ch, fil_bar_r_ch, bed_ch)
    umi_ch       = umi_combining(bam_r10_ch, bam_hd_ch, bed_ch, merged_ch, umi_r_for_ch, umi_r_rev_ch)
    bandu_ch     = final_filter(bar_fil_ch, umi_ch, bandu_r_for_ch, bandu_r_rev_ch)
    pairs_ch     = umi_pairing(bandu_ch, pairing_r_ch)
    degen_ch     = umi_degeneracy(bandu_ch, degen_r_ch)
    plot_degeneracy(degen_ch, plot_r_ch)
    summary_stats(pairs_ch, bandu_ch)
}

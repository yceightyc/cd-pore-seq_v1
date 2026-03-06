#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.bam_ribo   = "./input/ribo_filtered.bam"
params.bam_ctl    = "./input/control_filtered.bam"
params.ref        = "./input/reference.fasta"
params.btdesign   = "./input/barcode_tri_design.tsv"
params.roi_bed    = "./input/UMIs_randomer_barcode.bed"
params.joinbase_r = "./UMI_permut_combining_bases.R"
params.correct_py = "./correction.py"
params.plot_acc_r = "./plot_randomer_acc.R"
params.tri_design = "./input/contig_tridesign.tsv"
params.rsig_real  = "./input/eventalign_ribo/real_sig.tsv"
params.rsig_sim   = "./input/eventalign_ribo/sim_sig.tsv"
params.getinput_r = "./prepare_Minput.R"


process pileup {
    publishDir "data/"

    input:
    path bam_ribo
    path bam_ctl
    path ref
    path roi_bed

    output:
    path "mpileup_output_ribo.maf"
    path "mpileup_output_ctl.maf"

    script:
    """
    samtools mpileup $bam_ribo -f $ref --no-output-ins --no-output-del --no-output-ends \
      -l $roi_bed -o mpileup_output_ribo.maf --output-QNAME -A -B --max-depth 0
    samtools mpileup $bam_ctl -f $ref --no-output-ins --no-output-del --no-output-ends \
      -l $roi_bed -o mpileup_output_ctl.maf --output-QNAME -A -B --max-depth 0
    """
}

process split_bases {
    publishDir "data/base_information"

    input:
    path maf_ribo
    path maf_ctl

    output:
    path "baseATpositions_ribo", type: 'dir'
    path "baseATpositions_ctl", type: 'dir'

    script:
    """
    mkdir -p baseATpositions_ribo
    mkdir -p baseATpositions_ctl

    awk '{print "echo "\$5">baseATpositions_ribo/"\$1"_"\$2"_a"}' $maf_ribo | sh
    awk '{print "echo "\$5">baseATpositions_ctl/"\$1"_"\$2"_a"}' $maf_ctl | sh
    awk '{print "echo "\$1"\\t"\$2"\\t"\$7">baseATpositions_ribo/"\$1"_"\$2"_b"}' $maf_ribo | sh
    awk '{print "echo "\$1"\\t"\$2"\\t"\$7">baseATpositions_ctl/"\$1"_"\$2"_b"}' $maf_ctl | sh
    """
}

process mutate_and_merge {
    publishDir "data/base_information"

    input:
    path baseATpositions_ribo
    path baseATpositions_ctl

    output:
    path "baseATpositions_ribo/merged_mutated_ribo", type: 'dir'
    path "baseATpositions_ctl/merged_mutated_ctl", type: 'dir'

    script:
    """
    mkdir -p baseATpositions_ribo/a_mutated_ribo baseATpositions_ctl/a_mutated_ctl
    mkdir -p baseATpositions_ribo/b_mutated_ribo baseATpositions_ctl/b_mutated_ctl
    mkdir -p baseATpositions_ribo/merged_mutated_ribo baseATpositions_ctl/merged_mutated_ctl

    for file in baseATpositions_ribo/*a; do
        name=\$(basename "\$file")
        awk -v FS="" '{for (i=1;i<=NF;i++) printf \$i"\\n"}' "\$file" |
        awk '{if (\$0 != "+" && \$0 != "-" && \$0 !~ /^[0-9]*\$/) print \$0}' \
        > baseATpositions_ribo/a_mutated_ribo/\${name}_mutated
    done

    for file in baseATpositions_ribo/*b; do
        name=\$(basename "\$file")
        awk '{print \$3}' "\$file" | sed 's/,/\\n/g' \
        > baseATpositions_ribo/b_mutated_ribo/\${name}_mutated
    done

    for file in baseATpositions_ribo/*a; do
        id=\$(basename "\$file" | awk -F'_' '{print \$1"_"\$2}')
        paste baseATpositions_ribo/a_mutated_ribo/\${id}_a_mutated \
            baseATpositions_ribo/b_mutated_ribo/\${id}_b_mutated \
        > baseATpositions_ribo/merged_mutated_ribo/\${id}_merged
    done

    for file in baseATpositions_ctl/*a; do
        name=\$(basename "\$file")
        awk -v FS="" '{for (i=1;i<=NF;i++) printf \$i"\\n"}' "\$file" |
        awk '{if (\$0 != "+" && \$0 != "-" && \$0 !~ /^[0-9]*\$/) print \$0}' \
        > baseATpositions_ctl/a_mutated_ctl/\${name}_mutated
    done

    for file in baseATpositions_ctl/*b; do
        name=\$(basename "\$file")
        awk '{print \$3}' "\$file" | sed 's/,/\\n/g' \
        > baseATpositions_ctl/b_mutated_ctl/\${name}_mutated
    done

    for file in baseATpositions_ctl/*a; do
        id=\$(basename "\$file" | awk -F'_' '{print \$1"_"\$2}')
        paste baseATpositions_ctl/a_mutated_ctl/\${id}_a_mutated \
            baseATpositions_ctl/b_mutated_ctl/\${id}_b_mutated \
        > baseATpositions_ctl/merged_mutated_ctl/\${id}_merged
    done
    """
}

process merge_bases {
    publishDir "data/base_information"

    input:
    path roi_bed
    path joinbase_r
    path baseATpositions_ribo_merged
    path baseATpositions_ctl_merged

    output:
    path "barcodes_umis_randomers_baseinfo_ribo.tsv"
    path "barcodes_umis_randomers_baseinfo_ctl.tsv"

    script:
    """
    for r in \$(cut -f1 $roi_bed | sort -u); do
          Rscript $joinbase_r \$r $baseATpositions_ribo_merged
          Rscript $joinbase_r \$r $baseATpositions_ctl_merged
        done

    ls $baseATpositions_ribo_merged/*all.tsv|awk -F'[/|_]' '{print \$0"\t"\$4}'|while read a b;do awk '{print "'\$b'\t"\$0}' \$a;done|grep -wv rn|grep -v "*"|less> barcodes_umis_randomers_baseinfo_ribo.tsv ;
    ls $baseATpositions_ctl_merged/*all.tsv|awk -F'[/|_]' '{print \$0"\t"\$4}'|while read a b;do awk '{print "'\$b'\t"\$0}' \$a;done|grep -wv rn|grep -v "*"|less> barcodes_umis_randomers_baseinfo_ctl.tsv ;
    """
}

process randomer_base {
    publishDir "data/base_information"

    input:
    path baseATpositions_ribo_merged
    path baseATpositions_ctl_merged

    output:
    path "randomer_basecalling_accuracy_ribo", type: 'dir'
    path "randomer_basecalling_accuracy_ctl", type: 'dir'

    script:
    """
    mkdir -p randomer_basecalling_accuracy_ribo randomer_basecalling_accuracy_ctl ;
    echo 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79|sed 's/ /\\n/g'|awk '{print "cat $baseATpositions_ribo_merged/*_"\$0"*>randomer_basecalling_accuracy_ribo/"\$0"_for_baseofallreads"}'|sh 
    echo 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79|sed 's/ /\\n/g'|awk '{print "cat $baseATpositions_ctl_merged/*_"\$0"*>randomer_basecalling_accuracy_ctl/"\$0"_for_baseofallreads"}'|sh
    echo 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208|sed 's/ /\\n/g'|awk '{print "cat $baseATpositions_ribo_merged/*_"\$0"*>randomer_basecalling_accuracy_ribo/"\$0"_rev_baseofallreads"}'|sh 
    echo 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208|sed 's/ /\\n/g'|awk '{print "cat $baseATpositions_ctl_merged/*_"\$0"*>randomer_basecalling_accuracy_ctl/"\$0"_rev_baseofallreads"}'|sh 
    """
}

process getref {
    publishDir "data/", mode:"copy"

    input:
    path bt_design
    path base_info_ribo
    path base_info_ctl

    output:
    path "all_reads_ref_ribo.fasta"
    path "all_reads_ref_ctl.fasta"

    script:
    """
    awk 'BEGIN{FS="\t"} 
      NR==FNR {
      barcode_map[\$1]=\$2"\t"\$3"\t"\$4"\t"\$5;
      next}
      (\$1 in barcode_map) {
      split(barcode_map[\$1], arr, "\t")
      seq = "TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCTGCCGATCT" \
          arr[1] \
          "CTGATGCGTAGCT" \$6 \
          arr[2] \
          \$7 "GCGTGACATCGATGACAGATA" \
          \$8 "CTCGTATCAGCACACAATAGTCCATTATACGCGCGTATAATGGACTATTGTGTGCTGATACGAG" \
          \$3 "TATCTGTCATCGATGTCACGC" \
          \$4 arr[4] \
          \$5 "AGCTACGCATCAG" arr[3] \
          "AGATCGGCAGCAATACGTAACTGAACGAAGTACAGG";
      print ">"\$1"_"\$2"\\n"seq}' $bt_design $base_info_ribo > all_reads_ref_ribo.fasta
    
    awk 'BEGIN{FS="\t"} 
      NR==FNR {
      barcode_map[\$1]=\$2"\t"\$3"\t"\$4"\t"\$5;
      next}
      (\$1 in barcode_map) {
      split(barcode_map[\$1], arr, "\t")
      seq = "TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCTGCCGATCT" \
          arr[1] \
          "CTGATGCGTAGCT" \$6 \
          arr[2] \
          \$7 "GCGTGACATCGATGACAGATA" \
          \$8 "CTCGTATCAGCACACAATAGTCCATTATACGCGCGTATAATGGACTATTGTGTGCTGATACGAG" \
          \$3 "TATCTGTCATCGATGTCACGC" \
          \$4 arr[4] \
          \$5 "AGCTACGCATCAG" arr[3] \
          "AGATCGGCAGCAATACGTAACTGAACGAAGTACAGG";
      print ">"\$1"_"\$2"\\n"seq}' $bt_design $base_info_ctl > all_reads_ref_ctl.fasta
    """
}

process indexref {
    publishDir "data/", mode:"copy"

    input:
    path all_ref_ribo
    path all_ref_ctl

    output:
    path "all_reads_ref_ribo.fasta.fai"
    path "all_reads_ref_ctl.fasta.fai"

    script:
    """
    samtools faidx $all_ref_ribo
    samtools faidx $all_ref_ctl
    """
}

process indexbam {
    publishDir "input/", mode:"copy"

    input:
    path bam_ribo
    path bam_ctl

    output:
    path "ribo_filtered.bam.bai"
    path "control_filtered.bam.bai"

    script:
    """
    samtools index $bam_ribo
    samtools index $bam_ctl
    """
}

process correction {
    publishDir "data/"

    input:
    path fastas_ribo
    path fastas_ctl
    path bam_ribo
    path bam_ctl
    path correct_py

    output:
    path "corrected_reads_ribo", type: 'dir'
    path "corrected_reads_ctl", type: 'dir'

    script:
    """
    mkdir corrected_reads_ribo
    mkdir corrected_reads_ctl
    python $correct_py $fastas_ribo $bam_ribo corrected_reads_ribo
    python $correct_py $fastas_ctl $bam_ctl corrected_reads_ctl
    """
}

process pileup_correct {
    publishDir "data/"

    input:
    path corrected_ribo_dir
    path corrected_ctl_dir
    path fastas_ribo
    path fastas_ctl

    output:
    path "piled_up_ribo_corrected.maf"
    path "piled_up_ctl_corrected.maf"

    script:
    """
    mkdir pileup_ribo
    mkdir pileup_ctl
    samtools view $fastas_ribo|awk '{print \$1"\\t62\\t79\\n"\$1"\\t191\\t208"}' > randomer_r.bed 
    samtools view $fastas_ctl|awk '{print \$1"\\t62\\t79\\n"\$1"\\t191\\t208"}' > randomer_c.bed
    ls $corrected_ribo_dir|awk -F'.' '{print "samtools mpileup $corrected_ribo_dir/"\$0" -f $fastas_ribo --no-output-ins --no-output-del --no-output-ends -l randomer_r.bed -o pileup_ribo/"\$1".maf --output-QNAME -A -B --max-depth 0"}'| sh
    ls $corrected_ctl_dir|awk -F'.' '{print "samtools mpileup $corrected_ctl_dir/"\$0" -f $fastas_ctl --no-output-ins --no-output-del --no-output-ends -l randomer_c.bed -o pileup_ctl/"\$1".maf --output-QNAME -A -B --max-depth 0"}'| sh
    find pileup_ribo -type f -print0 | xargs -0 cat > piled_up_ribo_corrected.maf
    find pileup_ctl -type f -print0 | xargs -0 cat > piled_up_ctl_corrected.maf
    """
}

process plot_acc {
    publishDir "data/", mode: "copy"

    input:
    path plot_acc_r
    path bam_ribo
    path bam_ctl
    path raw_base_ribo_dir
    path raw_base_ctl_dir
    path corrected_maf_ribo
    path corrected_maf_ctl
    path contig_design

    output:
    path "randomer_acc_raw.pdf"
    path "randomer_acc_corrected.pdf"

    script:
    """
    samtools view $bam_ribo|cut -f1,3 > contig_info_ribo.tsv
    samtools view $bam_ctl|cut -f1,3 > contig_info_ctl.tsv
    Rscript $plot_acc_r $raw_base_ribo_dir $corrected_maf_ribo $contig_design $raw_base_ctl_dir $corrected_maf_ctl randomer_acc_raw.pdf randomer_acc_corrected.pdf
    """
}

process model_input {
    publishDir "data/"

    input:
    path dir_real_sig
    path dir_sim_sig
    path getinput_r

    output:
    path "model_input.tsv"

    script:
    """
    Rscript $getinput_r $dir_real_sig $dir_sim_sig model_input.tsv
    """
}

workflow {
    // Define channels for inputs
    bam_ribo_ch  = Channel.fromPath(params.bam_ribo)
    bam_ctl_ch   = Channel.fromPath(params.bam_ctl)   
    ref_ch       = Channel.fromPath(params.ref)
    roi_bed_ch   = Channel.fromPath(params.roi_bed)
    btdesign_ch  = Channel.fromPath(params.btdesign)
    join_r_ch    = Channel.fromPath(params.joinbase_r)
    corr_py_ch   = Channel.fromPath(params.correct_py)
    plot_acc_ch  = Channel.fromPath(params.plot_acc_r)
    tri_des_ch   = Channel.fromPath(params.tri_design)
    rsig_real_ch = Channel.fromPath(params.rsig_real)
    rsig_sim_ch  = Channel.fromPath(params.rsig_sim)
    getinputr_ch = Channel.fromPath(params.getinput_r)


    // Stepwise pipeline
    maf_ch        = pileup(bam_ribo_ch, bam_ctl_ch, ref_ch, roi_bed_ch)
    bases_ch      = split_bases(maf_ch) 
    merged_ch     = mutate_and_merge(bases_ch)
    base_join_ch  = merge_bases(roi_bed_ch, join_r_ch, merged_ch)
    ran_base_ch   = randomer_base(merged_ch)
    fastas_ch     = getref(btdesign_ch, base_join_ch)
    fastas_idx_ch = indexref(fastas_ch)
    bam_idx_ch    = indexbam(bam_ribo_ch, bam_ctl_ch)
    corrected_ch  = correction(fastas_ch, bam_ribo_ch, bam_ctl_ch, corr_py_ch)
    maf_cor_ch    = pileup_correct(corrected_ch, fastas_ch)
    plot_acc(plot_acc_ch, bam_ribo_ch, bam_ctl_ch, ran_base_ch, maf_cor_ch, tri_des_ch)
    m_input_ch    = model_input(rsig_real_ch, rsig_sim_ch, getinputr_ch)
}

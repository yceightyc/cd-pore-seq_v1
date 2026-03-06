#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.bams       = "./input/bams"
params.alnfastas  = "./input/alignment_fastas"
params.allref     = "./input/all.fasta"
params.bed        = "./input/UMIs_mods_oppo.bed"
params.cal_acc_r  = "./cal_duplex_acc.R"
params.eventalign = "./input/eventalign_all.tsv"
params.signal_r   = "./signal_analysis.R"
params.model_py   = "./model_training.py"

process bam2fastq {
    publishDir "data/"

    input:
    path bams

    output:
    path "fastqs", type: 'dir'

    script:
    """
    mkdir fastqs
    ls bams/|awk -F '.' '{print "samtools fastq -T \\"*\\" bams/"\$1".bam > fastqs/simplex_"\$1".fastq"}'|sh
    """
}

process minimap2_align {
    publishDir "data/"

    input:
    path fastqs
    path refs

    output:
    path "alignment", type: 'dir'

    script:
    """
    mkdir alignment
    ls $fastqs|awk -F '[_|.]' '{print "minimap2 --MD -y -a $refs/"\$3".fasta fastqs/"\$0" -k 5 -w 3|samtools sort -T tmp -o alignment/minimap2_"\$3".bam"}'|sh 
    """
}

process filter_primary {
    publishDir "data/"

    input:
    path bams

    output:
    path "filtered", type: 'dir'

    script:
    """
    mkdir filtered
    ls $bams/*.bam|awk -F'[_|.]' '{print "samtools view -bS -F 256 -F 16 -F 4 -F 2048 "\$0" | \
      samtools sort -T tMp -o filtered/"\$2"_forward.bam"}'|sh
    
    ls $bams/*.bam|awk -F'[_|.]' '{print "samtools view -bS -F 256 -f 16 -F 4 -F 2048 "\$0" | \
      samtools sort -T tmp -o filtered/"\$2"_reverse.bam"}'|sh
    """
}

process filter_concatemers {
    publishDir "data/filtered"

    input:
    path bams

    output:
    path "filtered/conca", type: 'dir'
 
    script:
    """
    mkdir -p filtered/conca
    cd filtered
    # The correct full size of a read should be 256 (including adapter seq) + 48 (barcode) = 304, set threshold as 320 to allow small insertions
    ls *.bam| \
      while read l;do samtools view \$l| \
      awk '{print \$1"\\t"\$4"\\t"length(\$10)}'| \
      awk '{print \$1"\\t"\$2+\$3}'| \
      awk '{if (\$2<=320) print \$1}' > "conca/\${l}_noconca_RN.rn"; done 
    """
}

process self_annealing {
    publishDir "data/filtered"

    input:
    path bams

    output:
    path "filtered/SA", type: 'dir'

    script:
    """
    mkdir -p filtered/SA
    cd filtered
    ls *.bam| \
      awk -F'.' '{print "samtools view "\$0"|grep \\"GCGTGACATCGATGACACGC\\"|cut -f1 > "SA/"\$1"_"\$2"_selfannealing_RN.rn"}'|sh    
    """
}

process intersect_reads {
    publishDir "data/filtered"

    input:
    path conca_dir
    path sa_dir

    output:
    path "intersect", type: 'dir'

    script:
    """
    mkdir intersect
    ls conca/*noconca_RN.rn| \
      awk -F '[/|_|.]' '{print "cat "\$0" SA/"\$2"_"\$3"_likely_selfanneal.rn SA/"\$2"_"\$3"_likely_selfanneal.rn | \
      sort|uniq -u > intersect/non_conca_no_selfannealing_RN_"\$2"_"\$3".rn"}'|sh
    """
}

process filter_bam {
    publishDir "data/"

    input:
    path bams
    path rns

    output:
    path "filtered_bams", type:'dir'
    
    script:
    """
    mkdir filtered_bams
    ls filtered/*.bam|awk -F '[/|_|.]' '{print "samtools view -N intersect/non_conca_no_selfannealing_RN_"\$2"_"\$3".rn -o filtered_bams/no_concat_noSA_"\$2"_"\$3".bam "\$0}'|sh
    """
}

process get_bed {
    publishDir "data/"

    input:
    path filtered_bams

    output:
    path "filtered_beds", type:'dir'

    script:
    """
    mkdir filtered_beds
    ls filtered_bams/*.bam|awk -F'[.|/]' '{print "bedtools bamtobed -i "\$0"|less > filtered_beds/"\$2".bed"}'|sh
    ls filtered_beds/no*reverse.bed|while read l;do awk '{print \$1"\t"198-\$3"\t"198-\$2"\t"\$4}' \$l > "\${l}_reversed";done
    ls filtered_beds/no*.bed*|awk -F '_' '{print \$5}'|sort -u|awk '{print "cat filtered_beds/no_concat_noSA_"\$0"_forward.bed filtered_beds/no_concat_noSA_"\$0"_reverse.bed_reversed|cut -f1-4 > filtered_beds/no_concat_noSA_"\$0"_revtransformed_all.bed"}'|sh ;
    """
}

process mod_truncation {
    publishDir "data/"

    input:
    path filtered_beds

    output:
    path "stats_truncation_atMods.tsv"

    script:
    """
    ls filtered_beds/no_concat_noSA_*_all.bed| \
      awk -F'_' '{print "less "\$0"|wc -l > filtered_beds/no_concat_noSA_"\$5"_readcount_all.txt"}'|sh
    ls filtered_beds/no_concat_noSA_*_all.bed| \
      while read l;do awk '{if (\$3>= 25 && \$3 <= 45) print \$0}' \$l|wc -l > "\${l}_stopnearmod.txt";done
    ls filtered_beds/no_concat_noSA_*_all.bed| \
      while read l;do awk '{if (\$3>=154 && \$3<=174) print \$0}' \$l|wc -l > "\${l}_stopopposemod.txt";done
    ls filtered_beds/no_concat_noSA_*all*.txt| \
      awk -F '_' '{print "paste filtered_beds/no_concat_noSA_"\$5"_readcount_all.txt filtered_beds/no_concat_noSA_"\$5"_revtransformed_all.bed_stopnearmod.txt filtered_beds/no_concat_noSA_"\$5"_revtransformed_all.bed_stopopposemod.txt > no_concat_noSA_"\$5"_modstopping_ratio.tsv"}'|sh
    ls no_concat_noSA_barcode*_modstopping_ratio.tsv|awk -F'_' '{print \$0"\\t"\$4}'|while read a b; do less filtered_beds/no_concat_noSA_"\${b}"_forward.bed|head -1|cut -f1|awk '{print "'\$a'\\t"\$0}' ;done|while read a b; do less \$a|awk '{print "'\$b'\\t"\$0"\t"\$2/\$1"\t"\$3/\$1}';done > stats_truncation_atMods.tsv
    """
}

process filter_trun {
    publishDir "data/", mode:'copy'

    input:
    path filtered_beds
    path bams

    output:
    path "allmods_filtered_final.bam"

    script:
    """
    ls filtered_beds/no_concat_noSA_barcode*_revtransformed_all.bed|while read l ; do less \$l|awk '{if (\$2 <= 26 && \$3 >= 173) print \$0}'|cut -f4 > "\${l}_essential.rn"; done 
    ls filtered_beds/*essential.rn|awk -F '_' '{print "samtools view -bS -F 256 -F 4 -F 2048 -N "\$0" -o filtered_final_all_"\$5".bam alignment/minimap2_"\$5".bam"}'|sh
    ls filtered_final_all_barcode* > all_filtered_bams.txt 
    samtools merge -b all_filtered_bams.txt -o allmods_filtered_final.bam
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
    publishDir "data/", mode:"move"

    input:
    path allref
    path merged_base
    path cal_UMI_match_r

    output:
    path "UMI_matched_rate_allmods.tsv"

    script:
    """
    mkdir UMI_permut_bases_out 
    
    samtools view $allref|cut -f1| \
      awk '{print "Rscript $cal_UMI_match_r "\$0" $merged_base"}'|sh
    
    ls UMI_permut_bases_out/*raw*all.tsv|awk -F'[/|_]' '{print \$0"\\t"\$6}'| \
      while read a b;do awk '{print "'\$b'\\t"\$0}' \$a;done| \
      grep -wv rn|grep -v "*" > umis_matching_aggregated.tsv
    
    # duplex rate per modification
    less umis_matching_aggregated.tsv|cut -f1|sort -u| \
      while read l;do awk '{if (\$1 == "'\$l'" && \$9 ==1 ) print \$0}' umis_matching_aggregated.tsv|wc -l| \
      awk '{print "'\$l'\\t"\$0}';done > UMI_matched_readcount_allmods.txt 
    less umis_matching_aggregated.tsv|cut -f1|sort -u| \
      while read l;do awk '{if (\$1 == "'\$l'") print \$0}' umis_matching_aggregated.tsv|wc -l| \
      awk '{print "'\$l'\\t"\$0}';done > UMI_all_readcount_allmods.txt
    paste UMI_matched_readcount_allmods.txt UMI_all_readcount_allmods.txt| \
      awk '{print \$1"\\t"\$2/\$4}' > UMI_matched_rate_allmods.tsv
    """
}

process mod_basecalls {
    publishDir "data/", mode:"move"

    input:
    path merged_base

    output:
    path "mod_basecalls_rate.tsv"
    path "oppomod_basecalls_rate.tsv"
    path "surrounding_mod_acc.tsv"
    path "surrounding_oppo_acc.tsv"

    script:
    """
    ls merged_mutated/*_35_merged | \
      while read l; do awk '{print toupper(\$1)}' \$l |sort|uniq -c|awk '{print "'\$l'\\t"\$0}';done| \
      awk -F '/' '{print \$2}'|awk -F'[\\t|_]' '{print \$1"\\t"\$4"\\t"\$5}' > mod_basecalls_readcount_allmods.txt 
    
    ls merged_mutated/*_164_merged | \
      while read l; do awk '{print toupper(\$1)}' \$l |sort|uniq -c|awk '{print "'\$l'\\t"\$0}';done| \
      awk -F '/' '{print \$2}'|awk -F'[\\t|_]' '{print \$1"\\t"\$4"\\t"\$5}' > oppomod_basecalls_readcount_allmods.txt
    
    awk '{sum[\$1]+=\$2; line[NR]=\$0} END {for(i=1;i<=NR;i++) {split(line[i],f);print line[i]"\\t"sum[f[1]]}}' mod_basecalls_readcount_allmods.txt| \
      awk '{print \$1"\\t"\$3"\\t"\$2/\$4}' > mod_basecalls_rate.tsv
    
    awk '{sum[\$1]+=\$2; line[NR]=\$0} END {for(i=1;i<=NR;i++) {split(line[i],f);print line[i]"\\t"sum[f[1]]}}' oppomod_basecalls_readcount_allmods.txt| \
      awk '{print \$1"\\t"\$3"\\t"\$2/\$4}' > oppomod_basecalls_rate.tsv

    sed -i 's/\\*/del/g' mod_basecalls_rate.tsv
    sed -i 's/\\*/del/g' oppomod_basecalls_rate.tsv

    # Quantify mod and opposite pos surrounding region basecalling accuracy
    for n in \$(seq -8 8); do m=\$((35+n)); for l in merged_mutated/*_\${m}_merged; 
      do [ -e "\$l" ] || continue; awk -v off="\$n" '{if(\$1=="."||\$1==",") sum++; count++} END {print FILENAME"\\t"sum"\\t"count"\\t"sum/count"\\t"count-sum"\\t"off}' "\$l"; done| \
      awk -F'/' '{print \$2}' | awk -F'[_|\\t]' '{print \$1"\\t"\$4"\\t"\$6"\\t'\$n'"}'; done > surrounding_mod_acc.tsv
    
    for n in \$(seq -8 8); do m=\$((164+n)); for l in merged_mutated/*_\${m}_merged; 
      do [ -e "\$l" ] || continue; awk -v off="\$n" '{if(\$1=="."||\$1==",") sum++; count++} END {print FILENAME"\\t"sum"\\t"count"\\t"sum/count"\\t"count-sum"\\t"off}' "\$l"; done| \
      awk -F'/' '{print \$2}' | awk -F'[_|\\t]' '{print \$1"\\t"\$4"\\t"\$6"\\t'\$n'"}' ; done > surrounding_oppo_acc.tsv
    """
}

process signal_analysis {
    publishDir "data/", mode:"copy"

    input:
    path signal_r
    path eventalign

    output:
    path "signal_analysis", type:"dir"

    script:
    """
    mkdir -p signal_analysis 
    
    Rscript $signal_r
    """
}

process model_train {
    publishDir "data/", mode:"copy"

    input:
    path model_py
    path model_input

    output:
    path "ROC_allmods.pdf"

    script:
    """
    python $model_py
    """
}

workflow {
    // Define channels for inputs
    bams_ch      = Channel.fromPath(params.bams, type: 'dir')
    alnref_ch    = Channel.fromPath(params.alnfastas, type: 'dir')  
    allref_ch    = Channel.fromPath(params.allref)
    roi_bed_ch   = Channel.fromPath(params.bed)
    cal_acc_r_ch = Channel.fromPath(params.cal_acc_r)
    event_ch     = Channel.fromPath(params.eventalign)
    signal_r_ch  = Channel.fromPath(params.signal_r)
    model_py_ch  = Channel.fromPath(params.model_py)

    // Stepwise pipeline
    fastq_ch     = bam2fastq(bams_ch)
    aln_ch       = minimap2_align(fastq_ch, alnref_ch)
    primary_ch   = filter_primary(aln_ch)
    conca_ch     = filter_concatemers(primary_ch)
    sa_ch        = self_annealing(primary_ch)
    intersect_ch = intersect_reads(conca_ch, sa_ch)
    bam_fil_ch   = filter_bam(primary_ch, intersect_ch)
    bed_ch       = get_bed(bam_fil_ch)
    mod_truncation(bed_ch)
    bam_final_ch = filter_trun(bed_ch, aln_ch)
    maf_ch       = pileup(bam_final_ch, allref_ch, roi_bed_ch)
    bases_ch     = split_bases(maf_ch)
    merged_ch    = mutate_and_merge(bases_ch)
    cal_pairing_acc(allref_ch, merged_ch, cal_acc_r_ch)
    mod_basecalls(merged_ch)
    signal_ch    = signal_analysis(signal_r_ch, event_ch)
    model_train(model_py_ch, signal_ch)
}
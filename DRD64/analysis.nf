#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.bam       = "./input/simplex.bam"
params.fastq     = "./input/simplex.fastq"
params.ref       = "./input/new_masked_ref.fasta"
params.unmaskref = "./input/unmask_reference_fasta_test10.fasta"
params.rawsignal = "./input/rawsignal/"
params.rscript   = "./plot_ampli_time_perkmer.R"
params.pyscript  = "./RF_training_val.py"
params.eventaligndir = "./"
params.sigdynar  = "./plot_ampli_dynamics.R"

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
    path "minimap2_masked_test.bam"

    script:
    """
    minimap2 --MD -y -a $ref $fastq -k 5 -w 3 | \
      samtools sort -T tmp -o minimap2_masked_test.bam
    """
}

process filter_primary_forward {
    publishDir "data/"

    input:
    path bam

    output:
    path "primary_for_minimap2_masked.bam"

    script:
    """
    samtools view -bS -F 256 -F 16 -F 4 -F 2048 $bam | \
      samtools sort -T tmp -o primary_for_minimap2_masked.bam
    """
}

process filter_concatemers {
    publishDir "data/"

    input:
    path bam

    output:
    path "noconca_RN_masked.rn"

    script:
    """
    samtools view $bam | \
      awk '{print \$1"\\t"\$4"\\t"length(\$10)}' | \
      awk '{print \$1"\\t"\$2+\$3}' | \
      awk '{if (\$2<=230) print \$0}' | \
      cut -f1 > noconca_RN_masked.rn
    """
}

process essential_regions {
    publishDir "data/"

    input:
    path bam

    output:
    path "essential_RN_masked.rn"

    script:
    """
    bedtools bamtobed -i $bam | \
      awk '{if (\$2 <= 55 && \$3 >= 110) print \$0}' | \
      cut -f4 > essential_RN_masked.rn
    """
}

process intersect_reads {
    publishDir "data/"

    input:
    path rn1
    path rn2

    output:
    path "non_conca_with_essential_RN_masked.rn"

    script:
    """
    cat $rn1 $rn2 | sort | uniq -d > non_conca_with_essential_RN_masked.rn
    """
}

process filter_bam {
    publishDir "data/"

    input:
    path bam
    path rn

    output:
    path "filterd_masked.bam"
    
    script:
    """
    samtools view -N $rn -o filterd_masked.bam $bam
    """
}

process index_filtered_bam {
    publishDir "data/"

    input:
    path bam_filtered

    output:
    path "filterd_masked.bam.bai"

    script:
    """
    samtools index $bam_filtered
    """
}


process f5c_eventalign {
    publishDir "data/"    

    input:
    path index
    path eventaligndir    

    output:
    path "../../../data/eventalign_output.tsv"

    script:
    """
    cd $eventaligndir ;
    sh eventalign.sh
    """
}

process plot_R {
    publishDir "data/"

    input:
    path rscript
    path eventalign
    path plotdir
	
    output:
    path "../../../data/model_input_final.tsv"

    script:
    """
    cd $plotdir ;
    Rscript $rscript
    """
}

process model_training {
    publishDir "data/"

    input:
    path pyscript
    path model_input
    path wrkdir
	
    output:
    path "../../../data/CM_plot.pdf"
    path "../../../data/PCA_plot.pdf"
    path "../../../data/multi-class_ROC.pdf"
    path "../../../data/Feature_Correlation.pdf"
    path "../../../data/Feature_Importance.pdf"

    script:
    """
    cd $wrkdir ;
    python $pyscript
    """
}

process plot_dynamics {
    publishDir "data/", mode:"move"

    input:
    path eventalign
    path dyna_r
	
    output:
    path "./signal_dyna/", type:'dir'

    script:
    """
    mkdir signal_dyna
    Rscript $dyna_r
    """
}

workflow {
    // Define channels for inputs
    bam_ch       = Channel.fromPath(params.bam)
    fastq_ch     = Channel.fromPath(params.fastq)
    ref_ch       = Channel.fromPath(params.ref)
    unmaskref_ch = Channel.fromPath(params.unmaskref)
    rawsignal_ch = Channel.fromPath(params.rawsignal)
    rscript_ch   = Channel.fromPath(params.rscript)
    pyscript_ch  = Channel.fromPath(params.pyscript)
    eventaligndir_ch     = Channel.fromPath(params.eventaligndir)  
    dyna_r_ch    = Channel.fromPath(params.sigdynar)

    // Stepwise pipeline
    fastq_ch     = bam2fastq(bam_ch)
    aln_ch       = minimap2_align(fastq_ch, ref_ch)
    primary_ch   = filter_primary_forward(aln_ch)
    conca_ch     = filter_concatemers(primary_ch)
    essential_ch = essential_regions(primary_ch)
    inter_ch     = intersect_reads(conca_ch, essential_ch)
    bam_filtered_ch      = filter_bam(primary_ch, inter_ch)
    index_ch     = index_filtered_bam(bam_filtered_ch)
    evt_ch       = f5c_eventalign(index_ch, eventaligndir_ch)
    rscript_ch   = plot_R(rscript_ch, evt_ch, eventaligndir_ch)
    model_training(pyscript_ch, rscript_ch, eventaligndir_ch)
    plot_dynamics(evt_ch, dyna_r_ch)
}

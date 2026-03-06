# Single DNA-embedded ribonucleotides Predictor
This repository contains a tool for **predicting single DNA-embedded ribonucleotides** from nanopore sequencing data.

**Input**:
``` 
  - Annotated raw signal (event-aligned output) saved in a folder; 
  - Reads alignment in .bam format
```

**Output**: 
```
  - Prediction result:
    Predicted nucleotide type samed in column 'rType_pred' (column 16), 0:dN; 1:rA; 2:rC; 3:rG; 4:rU
    Probability of being the predicted class in column 'rType_prob' (column 17)

```

---

## Quick Start
**Run Prediction:**
```bash
# Get forward and reverse names:
samtools view <alignment.bam>|awk '{if ($2 ==0) print $1"\tpos"}' > <primary_for.rn> ;

samtools view <alignment.bam>|awk '{if ($2 ==16) print $1"\trev"}'> <primary_rev.rn> ;

# Get bases at each position from reference sequence of alignment
samtools view <ref_seq.fasta> | cut -f1,10 | awk '{ref_pos=0; for(i=1;i<=length($2);i++){print $1,substr($2,i,1),ref_pos; ref_pos++}}' > <ref_base.tsv>;

# Run prediction
python predict_singleRibo.py <output_real_sig.tsv> <forward.rn> <reverse.rn> XGB_SRICO_1041.pkl <ref_base.tsv> <output_directory> ;

# Fetch predicted ribonucleotides and merge into single file
cat <output_directory>/*.txt|awk -F'\t' '{if ($16 != 0 && $1 != "read_name") print $0}' > <predicted_ribo.tsv>
```

i. e. 
```bash
samtools view pol2.bam|awk '{if ($2 ==0) print $1"\tpos"}' > pol2_forward.rn ;

samtools view pol2.bam|awk '{if ($2 ==0) print $1"\trev"}' > pol2_reverse.rn ;

samtools view reference.fasta | cut -f1,10 | awk '{ref_pos=0; for(i=1;i<=length($2);i++){print $1,substr($2,i,1),ref_pos; ref_pos++}}' > ref_base.tsv;

python predict_singleRibo.py event_aligned_pol2/ pol2_forward.rn pol2_reverse.rn XGB_SRICO_1041.pkl ref_base.tsv XGB_out_pol2

cat XGB_out_pol2/*.txt|awk -F'\t' '{if ($16 != 0 && $1 != "read_name") print $0}' > predicted_rN_pol2.tsv
```

---

## Table of Content
1. Example Input Annotated Signal Data (`./event_aligned_pol2/`)
2. Example XGB rNMP Prediction output (`./XGB_out_pol2/`)

---

## Required input
1. **Annotated Raw Signal** (generated with [f5c](https://github.com/hasindu2008/f5c))  
    - Formats: `.tsv`  
    - Example: `./event_aligned_pol2/*.txt`   

    Generate as follows, tools: [blue-crab](https://github.com/Psy-Fer/blue-crab); [SAMtools](https://github.com/samtools/samtools):
    ```bash
    # Convert pod5 into blow5 using blue-brab, recommended by f5c
    blue-crab p2s <pod5_directory> -o <real_signal.blow5>

    # Index basecalls.fastq and corrected_alignment.bam with signal.blow5 (GPU recommended)
    samtools fastq -T "*" <all_basecalls.bam> > <all_basecalls.fastq>
    f5c index --slow5 <real_signal.blow5> <all_basecalls.fastq>
    ## <all_basecalls.bam> could be dorado basecalling output, or filtered basecalls in .bam

    # Index alignment.bam by samtools
    samtools index <alignment.bam> 

    # Run f5c eventalignment (GPU recommended)
    f5c eventalign --slow5 <real_signal.blow5> --scale-events --print-read-names --pore r10 --min-recalib-events 20 -r <all_basecalls.fastq> -b <alignment.bam> -g <refseq.fasta> -o <output_real_sig.tsv>
    ## Output can be separated into chunks in the same folder to allow parallel processing`
    ```

2. **Forward read names** (from filtered alignment result, by [SAMtools](https://github.com/samtools/samtools))
    - Formats: `.rn`  
    - Example: `./pol2_forward.rn`

    Generate as follows:
    ```bash
    samtools view <alignment.bam>|awk '{if ($2 == 0) print $1"\tpos"}' > <forward.rn> 
    ```

3. **Reverse read names** (from filtered alignment result, by [SAMtools](https://github.com/samtools/samtools))
    - Formats: `.rn`  
    - Example: `./pol2_reverse.rn`

    Generate as follows:
    ```bash
    samtools view <alignment.bam>|awk '{if ($2 == 16) print $1"\trev"}' > <reverse.rn>
    ```

4. **Base at position on reference sequence** (from reference.fasta, by [SAMtools](https://github.com/samtools/samtools))
    - Formats: `.tsv`  
    - Example: `./ref_base.tsv` 
    
    (column 1: contig; column 2: base; column 3: position)

    Generate as follows:
    ```bash
    samtools view <ref_seq.fasta> | cut -f1,10 | awk '{ref_pos=0; for(i=1;i<=length($2);i++){print $1,substr($2,i,1),ref_pos; ref_pos++}}' > <ref_base.tsv>;
    ```

---

## To reproduce the yeast gDNA strand-enrichment analysis in the manuscript
(Replication origins are obtained from [OriDB](https://cerevisiae.oridb.org/))

run 
```bash
# Run prediction as above, then
python strand_analysis.py
```


---

## Reference
Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021).
Twelve years of SAMtools and BCFtools.
GigaScience, 10(2), giab008.
https://doi.org/10.1093/gigascience/giab008

Gamaarachchi, H., Lam, C. W., Jayatilaka, G., Samarakoon, H., Simpson, J. T., Smith, M. A., & Parameswaran, S. (2020).
GPU accelerated adaptive banded event alignment for rapid comparative nanopore signal analysis.
BMC Bioinformatics, 21(1), 343.
https://doi.org/10.1186/s12859-020-03659-2

Gamaarachchi, H., Samarakoon, H., Jenner, S. P., Ferguson, J. M., Amos, T. G., Hammond, J. M., Saadat, H., Smith, M. A., Parameswaran, S., & Deveson, I. W. (2022).
Fast nanopore sequencing data analysis with SLOW5.
Nature Biotechnology, 1–4.
https://doi.org/10.1038/s41587-022-01383-6

Reijns, M. A. M., Kemp, H., Ding, J., de Procé, S. M., Jackson, A. P., & Taylor, M. S. (2015). 
Lagging-strand replication shapes the mutational landscape of the genome. 
Nature, 518(7540), 502–506. 
https://doi.org/10.1038/nature14183

Siow, C. C., Nieduszynska, S. R., Müller, C. A., & Nieduszynski, C. A. (2012). 
OriDB, the DNA replication origin database updated and extended. 
Nucleic acids research, 40(D1), D682-D686.
https://doi.org/10.1093/nar/gkr1091
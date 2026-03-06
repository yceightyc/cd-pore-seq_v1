# CD-Pore-Seq_v1

CD-pore-seq is an ONT nanopore-based sequencing technique for generating context-informed high-complexity duplex data for DNA modification models; and offers a tool to detect single DNA-embedded ribonucleotides.

This repository contains:
1. **Scripts** and **input data** to reproduce the results in the manuscript
2. A **Tool with instructions** to predict **single DNA-embedded ribonucleotides** from nanopore sequencing

---
## Table of Content
Reproducing the results:
1. Fig. 1 & 2; Extended Data Fig. 1-3; Extended Data Table 1 (`./DRD64/`)
2. Extended Data Fig. 4 (`./ArUC/`)
3. Fig. 3b (`./UMI/`)
4. Fig. 3d & e; Extended Data Fig. 9 (`./TelN/`)
5. Fig. 3f; Fig. 4; Extended Data Fig. 5 & 6 (`./model_training/`)
6. Fig. 5b-d; Fig. 6; Extended Data Fig. 7 & 8; Extended Data Table 3 (`./Mixed_Modifications/`)
7. Single DNA-embedded ribonucleotides predictor (`./Single_rNMP_predictor/`); 

## Quick Start

- Build enviornment using [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
```bash
# Clone this repository
git clone https://git.ecdf.ed.ac.uk/s2404760/cd-pore-seq_v1.git 

# Build enviornment
cd cd-pore-seq_v1
conda env create -f conda.yml
```

**Download and extract input data from Edinburgh Datashare:**
```bash
# Go to the specific folder:
e.g.
cd DRD64/input

# Download respective data then unzip:
e.g.
tar -xvzf DRD64.tar.gz

# Then follow instructions to run related analysis.
```

**Predicing single DNA-embedded ribonucleotides from Nanopore (R10.4.1) sequencing:**
```bash
# GPU required
cd cd-pore-seq_v1
cd Single_rNMP_predictor
# Following respective instructions in the folder "Single_rNMP_predictor".
```
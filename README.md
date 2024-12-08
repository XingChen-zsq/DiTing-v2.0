# Metasnake
## Introduction
An easy-to-use wrapper for a robust snakemake pipeline designed for metagenomic assembly, gene prediction, genomic binning, taxonomy and functional annotations analysis of biogeochemical cycles. The pipeline features a series of flexible, interdependent modules optimized for fast execution. Each module can also be executed independently as a single pipeline, depending on the provided input. Interestingly, Diting 2.0v can be able to modify and add new modules based on snakemake workflow.
## Workflow
This workflow defines data analysis in terms of rules that are specified in the Snakefile.
![workflow](https://github.com/XingChen-zsq/Metasnake/blob/main/workflow.png)
## How to use
This workflow provides a series of .yaml files (including config.yaml and dependent environment configuration file, like kegg.yaml), module files () and one Snakefile. You can run  by editing the parameter file config.yaml  and the main file Snakefile.

 - **config.yaml:** including defines of input and output `folders`, `thread` and `cpu` configuration for specific rules. Before starting, you can specify or define the `input and output folder`, the number of `threads` and the `cpu` for specific running rules.
     ```
    # Input and output folder
     READS_SUF: ".fastq"
    READS_DIR: "cleanreads"
    ASSEMBLY_DIR: "Assembly"
    PRODIGAL_DIR: "Prodigal"  
    CDHIT_DIR: "CD_hit"  
    FILTERED_DIR: "Filtered"
    BBMAP_DIR: "BBMap"  
    BWA_INDEX: "BBMap/bwa_index"
    MAPPING: "BBMap/mapping"
    PILEUP: "BBMap/pileup"
    GENE_ABUN_DIR: "Abundance" 
    KODB_DIR: "kofam_database"        
    KEGG_DIR: "KEGG_annotation"       
    OUT_DIR: "final_output"
    TABLE: "table" 
    # Thread configuration for specific rules
    threads:
       megahit: 8
       prodigal: 4
       cdhit: 4
       bwa: 8
       metawrap_binning: 8
       metawrap_refinement: 8
       quantify_bins: 8
   # CPU configuration for specific rules
   cpu:
       kegg_annotation: 4
       gtdb_classification: 4
 - **Snakefile:** including `rule all` and a series of `modules`.
```
# Define the main workflow
rule all:
    input:
        expand(os.path.join(config["GENE_ABUN_DIR"], "{basename}.abundance"), basename=config["BASENAMES"]),
        os.path.join(config["KEGG_DIR"], 'pathways_relative_abundance_gene_level.tab'),
        os.path.join(config["OUT_DIR"], 'carbon_cycle_sketch.png'),
        os.path.join(config["OUT_DIR"], 'carbon_cycle_heatmap.pdf'),
        os.path.join(config["BIN_ABUNDANCE_DIR"], "bin_abundance_table.tab")

# input modules
module Assembly_gene_prediction:
    snakefile: "./Assembly_gene_prediction.smk"
    config: config

module Gene_abundance:
    snakefile: "./Gene_abundance.smk"
    config: config

module Function_annotation:
    snakefile: "./Function_annotation.smk"
    config: config

module Pathway_abundance:
    snakefile: "./Pathway_abundance.smk"
    config: config

module Visualization:
    snakefile: "./Visualization.smk"
    config: config

module Binning:
    snakefile: "./Binning.smk"
    config: config
    
# execute modules
use rule * from Assembly_gene_prediction as *
use rule * from Gene_abundance as *
use rule * from Function_annotation as *
use rule * from Pathway_abundance as *
use rule * from Visualization as *
use rule * from Binning as *
```
## Running

 

 **1. Download databases**
```
# At the home directory of this program
  mkdir kofam_database/
  cd kofam_database/
  wget -c ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz 
  wget -c ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz 
  gzip -d ko_list.gz
  tar zxvf profiles.tar.gz
```
 **2. Create input folder contain your data**
```
# At the home directory of this program
  mkdir cleanreads/ 
```
**3. Activate snakemake environment**
```
  conda activate snakemake
``` 

**4. After defining the `threads` or `cpu` of specific rules at `config.yaml` file, and choosing the rule all in`modules` of the specified rule, then run:**
```
  snakemake --use-conda --core 8
```
## Dependencies

 - **metagenome.yaml:** configuration environment of metagenomic analysis.
```
name: metagenome_env
channels:
 - bioconda
 - conda-forge
dependencies:
 - megahit >= 1.1.3
 - Prodigal >= 2.6.3
 - cd-hit >= 4.8.1
 - bwa >= 0.7.17
 - python>=3.8
 - hmmer
 - parallel
 - kofamscan
```
 - **genomic.yaml:** configuration environment of genomic binning.
```
name: genomic_env
channels:
 - ursky
 - bioconda
 - conda-forge
 - defaults
dependencies:
 - python=2.7
 - bwa
 - samtools
 - metabat2
 - concoct
 - MaxBin2
 - checkm-genome
 - prodigal
 - pplacer
 - numpy
 - matplotlib
 - pysam
 - hmmer
 - SPAdes
 - salmon
 - seaborn
 - metawrap-mg
```
## Download databases
 - **kofam_database**
 -   ko_list.gz (ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz)
-   profiles.tar.gz (ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz)

## Installation
      
 

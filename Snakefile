import os
import csv
import logging
from glob import glob

# Basic logging configuration
logging.basicConfig(level=logging.INFO)

# Input configuration files
configfile: "config.yaml"

# Get all sample base names
READS_FILES_1 = glob(os.path.join(config["READS_DIR"], f"*_1{config['READS_SUF']}"))
READS_FILES_2 = glob(os.path.join(config["READS_DIR"], f"*_2{config['READS_SUF']}"))

# Extract sample basenames
BASENAMES = [os.path.basename(f).replace(f'_1{config["READS_SUF"]}', '') for f in READS_FILES_1]

# input BASENAMES to config.yaml
config["BASENAMES"] = BASENAMES

# Define the main workflow
rule all:
    input:
        expand(os.path.join(config["GENE_ABUN_DIR"], "{basename}.abundance"), basename=config["BASENAMES"]),
        os.path.join(config["KEGG_DIR"], 'pathways_relative_abundance_gene_level.tab'),
        os.path.join(config["OUT_DIR"], 'carbon_cycle_sketch.png'),
        os.path.join(config["OUT_DIR"], 'carbon_cycle_heatmap.pdf'),
        os.path.join(config["BIN_ABUNDANCE_DIR"], "bin_abundance_table.tab")

## input module
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

use rule * from Assembly_gene_prediction as *
use rule * from Gene_abundance as *
use rule * from Function_annotation as *
use rule * from Pathway_abundance as *
use rule * from Visualization as *
use rule * from Binning as *



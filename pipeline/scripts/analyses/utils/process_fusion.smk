## Main snakemake for genewise analysis on genomic data
## Date : Decembre 2024
## Authors :


# Import python modules
# ---------------------
import os
import json

# Function to load JSON files
# ---------------------------
def load_json(file_path):
    with open(file_path, "r") as file:
        return json.load(file)

# Assign environment variables
# ----------------------------
globals().update(load_json("../environment_path.json"))


# Configuration
# -------------
configfile: "analyse.json"
configfile: "assemblies.json"


# List of assemblies
# ------------------
ACCESSNB = config["assembly_list"]

# List of domains to be fully processed. 
# --------------------------------------
DOMAINS = config["domains"]

# List of domains to be processed without paralogy checks. 
# --------------------------------------------------------
DOMAINS_SIMPLE = config["domains_simple"]

ALL_DOMAINS = DOMAINS + DOMAINS_SIMPLE

# The reference alignment of each domain 
# ---------------------------------------
DOMAIN_REFERENCES = config["domain_references"]

EXONS = config["exons"]

# Name of the resources directory. The name is  used  pathGTDriftResource:
# it contains the resources used by the analyse 
# The directories are structured as follows
# RESOURCES_DIR_NAME/

PROTEIN_RESOURCES_DIR_NAME  = config["protein_resources_dir_name"]

GENEWISE_RESOURCES_DIR_NAME  = config["genewise_resources_dir_name"]

# Name of global results directory.
# The directory is located in pathGTDriftGlobalResults
# ----------------------------------------------------
GLOBAL_RESULTS = config["analyse_dir_name"]

# Name of genome specific results directory.
# The directory is located in genome_assembly/{accession}/analyses/
# -----------------------------------------------------------------
GENOME_RESULTS = config["analyse_dir_name"]


# function to get the name of the reference alignment for a domain
# ----------------------------------------------------------------
def get_domain(wildcards):
    domain = wildcards.domain
    return domain 
    
# function to get the name of the reference alignment for a domain
# ----------------------------------------------------------------
def get_reference(wildcards):
    fname = DOMAIN_REFERENCES.get(wildcards.domain, "")
    return fname 
    
# function to get the path of the reference alignment for a domain
# ----------------------------------------------------------------
def get_reference_file(wildcards):
    fname = DOMAIN_REFERENCES.get(wildcards.domain, "")
    domain = wildcards.domain
    return pathGTDriftResource + RESOURCES_DIR_NAME + "hmm_profiles/"+ domain+"/"+fname+".hmm"
        
# -----------------------------------------------
# all : inputs define the to be files generated . 
# -----------------------------------------------       
            
rule all:
    input:
        protein_with_chromo=expand(pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS
        +"whole_summary_with_chromosomes.csv",accession=ACCESSNB), 
        protein_and_genomic=expand(pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS
        +"whole_summary_protein_genomic.csv",accession=ACCESSNB), 
        
rule get_proteins_chromosome_info:
    input:
        summary = pathGTDriftData
        + "genome_assembly/{accession}/analyses/"+PROTEIN_RESOURCES_DIR_NAME + "whole_summary.csv",
        protein = pathGTDriftData + "genome_assembly/{accession}/annotation/protein.faa", 
        gff = pathGTDriftData + "genome_assembly/{accession}/annotation/genomic.gff", 
    output:
        pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +"whole_summary_with_chromosomes.csv",
    shell:
        """
        python3 ../utils/python/get_chromosome_positions.py -i {input.summary}  -p {input.protein} -g {input.gff} -o {output}
        """
        
rule fusion_protein_genomic:
    input:
        genomic = pathGTDriftData
        + "genome_assembly/{accession}/analyses/"+GENEWISE_RESOURCES_DIR_NAME + "parsed_whole_summary_genewise.csv",
        protein = pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +"whole_summary_with_chromosomes.csv",
    output:
        pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +"whole_summary_protein_genomic.csv",
    shell:
        """
        python3 ../utils/python/fusion_protein_genewise.py -g {input.genomic}  -p {input.protein}  -o {output}
        """
                     
# Modules snakemake
# -----------------

#include: "../utils/module_genewise.smk"
#include: "../utils/module_check_paralogs_genewise.smk"
#include: "../utils/module_concatenate.smk"

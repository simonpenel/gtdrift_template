## snakemake for analysis the dna sequences of the zinc fingers protein detected with process_genewise.smk

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




RESOURCES_DIR_NAME  = config["resources_dir_name"]

# The reference alignment of each domain
# ---------------------------------------
DOMAIN_REFERENCES = config["domain_references"]

# List of domains to be fully processed.
# --------------------------------------
DOMAINS = config["domains"]

# List of domains to be processed without paralogy checks.
# --------------------------------------------------------
DOMAINS_SIMPLE = config["domains_simple"]

# Name of global results directory.
# The directory is located in pathGTDriftGlobalResults
# ----------------------------------------------------
GLOBAL_RESULTS = config["analyse_dir_name"]

# Name of genome specific results directory.
# The directory is located in genome_assembly/{accession}/analyses/
# -----------------------------------------------------------------
GENOME_RESULTS = config["analyse_dir_name"]

storagetype  = config["storagetype"]




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

# get the files and directories describing the reference alignments
# -----------------------------------------------------------------
directories, files = glob_wildcards(pathGTDriftResource + RESOURCES_DIR_NAME + "reference_alignments/{dir}/{file}.fst")

# Rules
# -----

# -----------------------------------------------
# all : inputs define the to be files generated .
# -----------------------------------------------
rule all:
    input:
        zincfinger_motif = expand(pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +  "zinc_finger_dna/zf_results.csv",accession=ACCESSNB),
        fasta_list = expand(pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +  "zinc_finger_dna/clustering/list_of_files.txt",accession=ACCESSNB),
 


# -------------------------------------------------------
# zincfinger_motif
# analyse of the genewise ouput
# -------------------------------------------------------
rule zincfinger_motif:
    """
    Run the zinc finger dna analysis on genewise output.
    """
    input:
        # Genewise output
        # ----------------
        genewise = "results/{accession}/Step3_genewise/gw.concat",
    output:
        results = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +  "zinc_finger_dna/zf_results.csv",
        log = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +  "zinc_finger_dna/zf_results.log",
        fasta = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +  "zinc_finger_dna/zf_results.fasta",
    shell:
        """
        python3 ../utils/python/genewise_parser_zf.py  -i {input.genewise} -a {wildcards.accession} -o {output.results} -f {output.fasta} -l {output.log}
        """      
       #python3 ../utils/python/get_ZF_motif_positions_in_protein.py -i {input.protein_seq} -g {input.gff}  -d {input.fasta_dna} -o {output.results} -l {output.log} -w {output.warnings}


# -------------------------------------------------------
# zincfinger_motif
# analyse of the dna sequences of the zincfinger proteins
# -------------------------------------------------------
rule get_fasta:
    """
    get the fasta file of zinc fingers.
    """
    input:
        results = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +  "zinc_finger_dna/zf_results.csv",
    output:
        fasta_list = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +  "zinc_finger_dna/clustering/list_of_files.txt",
    shell:
        """
        python3 ../utils/python/extract_zf_fasta.py -i {input.results} -o {pathGTDriftData}genome_assembly/{wildcards.accession}/analyses/{GENOME_RESULTS}zinc_finger_dna/clustering
        """
        
        
        
        

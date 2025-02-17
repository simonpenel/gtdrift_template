## Main snakemake for domain analysis on protein data
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
        zincfinger_out = expand(pathGTDriftGlobalResults + GLOBAL_RESULTS + "zinc_finger/zinc_finger_details/ZFD_{accession}.csv",accession=ACCESSNB),
        zinc_finger_combined = pathGTDriftGlobalResults + GLOBAL_RESULTS + "zinc_finger/zinc_finger.csv",
        global_zf_count=pathGTDriftGlobalResults 
        + GLOBAL_RESULTS + "zinc_finger/zf_count.csv",

rule zincfinger_analysis:
    """
    Run the zinc finger analysis on each protein sequence using R.
    """
    input:
        # Candidates after paralog checks
        # --------------------------------
        protein_seq = pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS
            + "selected_candidates.fasta",
        #protein_seq = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/candidates_prdm9.fasta"
    output:
        zincfinger_out = pathGTDriftGlobalResults + GLOBAL_RESULTS + "zinc_finger/zinc_finger_details/ZFD_{accession}.csv"
    
    run:
        # Verify if the input file exists
        if not os.path.exists(input.protein_seq):
            raise FileNotFoundError(f"Input file does not exist: {input.protein_seq}")
        
        # Run the R script
        command = f"Rscript --vanilla ../utils/zincfinger_analysis.R {input.protein_seq} {output.zincfinger_out}"
        shell(command)

rule combine_zinc_finger:
    """
    Combine all ZFD_{accession}.csv files into one zinc_finger.csv.
    """
    input:
        zfd_files = expand(pathGTDriftGlobalResults + GLOBAL_RESULTS + "zinc_finger/zinc_finger_details/ZFD_{accession}.csv", accession=ACCESSNB)
    output:
        zinc_finger_combined = pathGTDriftGlobalResults + GLOBAL_RESULTS + "zinc_finger/zinc_finger.csv"
    run:
        # Check if there are input files before trying to combine
        if not input.zfd_files:
            raise FileNotFoundError("No input files found for combining.")
        
        # Combine all ZFD_{accession}.csv files into one zinc_finger.csv
        with open(output.zinc_finger_combined, 'w') as outfile:
            for fname in input.zfd_files:
                with open(fname) as infile:
                    outfile.write(infile.read())
                    outfile.write("\n")  # Add a newline between files


# -----------------------------------------------------------------
# create_global_zf_table
# write a summary of zf protein with more than 5 zf domains.
# Output format:
# Accession;5+ ZF
# -----------------------------------------------------------------
rule create_global_zf_table:
    """
    Creation of global zf count table
    """
    input:
        expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/" 
            + GENOME_RESULTS 
            + "hmm_search/domtbl/ZF_domains_summary",
            accession=ACCESSNB,
        ),
    output:
        pathGTDriftGlobalResults 
        + GLOBAL_RESULTS + "zinc_finger/zf_count.csv",
    params:
        path=pathGTDriftData + "genome_assembly"      
    script:
         "../utils/python/zf_analysis.py"    
                     

# Modules snakemake
# -----------------



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

RESOURCES_DIR_NAME  = config["resources_dir_name"]


# Name of global results directory.
# The directory is located in pathGTDriftGlobalResults
# ----------------------------------------------------
GLOBAL_RESULTS = config["analyse_dir_name"]

# Name of genome specific results directory.
# The directory is located in genome_assembly/{accession}/analyses/
# -----------------------------------------------------------------
GENOME_RESULTS = config["analyse_dir_name"]

# config for hmm
# --------------
if config["mode"] == "guix":
    RUNCMD = "guix shell hmmer -- "
else:
    RUNCMD = ""
    
storagetype  = config["storagetype"]

if config["Nb_aln_genewise"] == "":
    ALN = '10'
else:
    ALN = config["Nb_aln_genewise"]

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
    print("wildcard "+str(wildcards))
    fname = DOMAIN_REFERENCES.get(wildcards.domain, "")
    domain = wildcards.domain
    print("domain " + domain + ", fname "+fname)
    return pathGTDriftResource + RESOURCES_DIR_NAME + "hmm_profiles/"+ domain+"/"+fname+".hmm"
        
# -----------------------------------------------
# all : inputs define the to be files generated . 
# -----------------------------------------------
rule all:
    """
    """
    input:
        #debug1=expand(pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/tbl/SET",accession=ACCESSNB),
        #debug2=expand(pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/domtbl/SET_domains",accession=ACCESSNB),
        debug1=expand(pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/tbl/{domain}",accession=ACCESSNB,domain=ALL_DOMAINS),
        debug2=expand(pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/domtbl/{domain}_domains",accession=ACCESSNB,domain=ALL_DOMAINS),

        # Statistics on candidates with all domains for all genomes
        # ---------------------------------------------------------           
        ##all_candidates_stats_summary=pathGTDriftGlobalResults + GLOBAL_RESULTS + "candidate_statistics_summary.csv",
        #sum=expand(pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +"whole_summary_genewise.csv",accession=ACCESSNB),  
        
        
        #parsed = expand(pathGTDriftData
        #+ "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +"parsed_whole_summary_genewise.csv",accession=ACCESSNB), 
        
        #wad = expand(pathGTDriftData
        #+ "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +"wad_parsed_whole_summary_genewise.csv",accession=ACCESSNB), 
        
        # Concatenation of results on all genomes
        #concat_parsed=pathGTDriftGlobalResults + GLOBAL_RESULTS + "ordered_parsed_results.csv",   
 
        # Concatenation of results on all genomes
        #concat_wad=pathGTDriftGlobalResults + GLOBAL_RESULTS + "wad_parsed_results.csv",           
        
        ## For the zincfinger analysis
        # candidates_confirmed=expand(pathGTDriftData 
        # + "genome_assembly/{accession}/analyses/" 
        # + GENOME_RESULTS + "candidates_SET.csv",accession=ACCESSNB), 
         
        # candidates_simple=expand(pathGTDriftData 
        # + "genome_assembly/{accession}/analyses/" 
        # + GENOME_RESULTS + "candidates_simple_ZF.txt",accession=ACCESSNB),           
         
        ## For the SET Tyronsine analysis
        # SET_fasta=expand(pathGTDriftData 
        # + "genome_assembly/{accession}/analyses/" 
        # + GENOME_RESULTS + "candidates_SET.fasta",accession=ACCESSNB), 
                    
# Modules snakemake
# -----------------

include: "../utils/module_genewise.smk"
include: "../utils/module_check_paralogs_genewise.smk"
include: "../utils/module_concatenate.smk"

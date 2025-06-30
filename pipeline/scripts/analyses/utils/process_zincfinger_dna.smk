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
configfile: "analyse_zf.json"
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
        zincfinger_hmm = expand(pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +  "zinc_finger_dna/ZF_hmm_{accession}.csv",accession=ACCESSNB),
        zincfinger_hmm_dna = expand(pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +  "zinc_finger_dna/ZF_hmm_{accession}_dna_id.csv",accession=ACCESSNB),
        zincfinger_analyse_out =  expand(pathGTDriftGlobalResults + GLOBAL_RESULTS + "zinc_finger_dna/zinc_finger_details/ZFD_{accession}.csv",accession=ACCESSNB),
        zincfinger_motif = expand(pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +  "zinc_finger_dna/zf_results.csv",accession=ACCESSNB),       


rule get_genome_seq_fasta:
    input:
        fasta = pathGTDriftData+ "genome_assembly/{accession}/genome_seq/genomic.fna.path"
    output:
        fasta = temp(pathGTDriftData+ "genome_assembly/{accession}/genome_seq/genomic.fna")
    shell:
        """
        export  genomic=`cat {input.fasta}`
        echo "Genome sequence fasta file : $genomic"
        if [ {storagetype} == irods ];
            then
            echo "iget  /lbbeZone/home/penel/gtdrift/genome_seq/$genomic"
            ls {pathGTDriftData}"genome_assembly/{wildcards.accession}/genome_seq/"
            iget -f /lbbeZone/home/penel/gtdrift/genome_seq/$genomic {output.fasta}
        else
            ln -s {pathGTDriftData}"genome_assembly/{wildcards.accession}/genome_seq/$genomic {output.fasta}"
        fi    
        """
        
# ------------
# get_hmm_zf :
# ------------
rule get_hmm_zf:
    input:
          ZF_per_domain=pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/domtbl/ZF_domains",                  
    output:
        zincfinger_hmm = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +  "zinc_finger_dna/ZF_hmm_{accession}.csv",
    shell:
        """
        cat {input.ZF_per_domain}| awk '{{print $1";"$18";"$19}}' | grep -v "#" > {output.zincfinger_hmm}
        """
# ----------------
# get_hmm_zf_dna :
# ----------------
rule get_hmm_zf_dna:
    input:
        zincfinger_hmm = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +  "zinc_finger_dna/ZF_hmm_{accession}.csv",  
        gff = pathGTDriftData + "genome_assembly/{accession}/annotation/genomic.gff",                
    output:
        zincfinger_hmm_dna_id = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +  "zinc_finger_dna/ZF_hmm_{accession}_dna_id.csv",
    shell:
        """
        python3 ../utils/python/get_ZF_positions_in_chromosome.py -i {input.zincfinger_hmm} -g {input.gff} -o {output}
        """

rule zincfinger_analysis:
    """
    Run the zinc finger analysis on each protein sequence using R.
    """
    input:
        # Candidates after paralog checks
        # --------------------------------
        protein_seq = pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS
            + "candidates_for_zf_analysis.fasta",
        #protein_seq = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/candidates_prdm9.fasta"
    output:
        zincfinger_out = pathGTDriftGlobalResults + GLOBAL_RESULTS + "zinc_finger_dna/zinc_finger_details/ZFD_{accession}.csv"
    
    run:
        # Verify if the input file exists
        if not os.path.exists(input.protein_seq):
            raise FileNotFoundError(f"Input file does not exist: {input.protein_seq}")
        
        # Run the R script
        command = f"Rscript --vanilla ../utils/zincfinger_analysis_dna.R {input.protein_seq} {output.zincfinger_out}"
        shell(command)

rule zincfinger_motif:
    """
    Run the zinc finger analysis on each protein sequence.
    """
    input:
        # Candidates after paralog checks
        # --------------------------------
        protein_seq = pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS
            + "candidates_for_zf_analysis.fasta",
        gff = pathGTDriftData + "genome_assembly/{accession}/annotation/genomic.gff", 
        fasta_dna = pathGTDriftData+ "genome_assembly/{accession}/genome_seq/genomic.fna",
        #protein_seq = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/candidates_prdm9.fasta"
    output:
        results = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +  "zinc_finger_dna/zf_results.csv",
        log = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +  "zinc_finger_dna/zf_results.log",
        warnings = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +  "zinc_finger_dna/zf_results.warnings",            
    shell:
        """
        python3 ../utils/python/get_ZF_motif_positions_in_protein.py -i {input.protein_seq} -g {input.gff}  -d {input.fasta_dna} -o {output.results} -l {output.log} -w {output.warnings}
        """
      
#    script:
#        "../utils/python/lol.py" 
        




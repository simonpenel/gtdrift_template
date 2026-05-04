## Module for prdm9 analysis on protein data
## Date : Decembre 2024
## Authors :


# Configuration
# -------------


configfile: "config.json"


ACCESSNB = config["assembly_list"]
DOMAIN = ["KRAB", "SET", "SSXRD", "ZF"]
#GLOBAL_RESULTS = "prm9_protein_analysis/"

# Function to load JSON files
# ---------------------------
def load_json(file_path):
    with open(file_path, "r") as file:
        return json.load(file)


# Assign environment variables
# ----------------------------
globals().update(load_json("../environment_path.json"))



# Rules
# -----

# -----------------------------------------------------------------
# create_global_krab_table
# write a summary of krab domain in the proteome.
# Output format:
# Accession;Protein List;KRAB nb
# -----------------------------------------------------------------
rule create_global_krab_table:
    """
    Creation of global krab result table 
    """
    input: 
        expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/tbl/KRAB_tabulated",
            accession=ACCESSNB,
        ),
    output:
        pathGTDriftGlobalResults + GLOBAL_RESULTS + "table_results/krab_data.csv",
    params:
        path=pathGTDriftData + "genome_assembly"
    
    script:
         "python/krab.py"
        
# -----------------------------------------------------------------
# create_global_krabzf_table
# write a summary of krab domain in the proteome.
# Output format:
# Accession;KRAB+ZF protein list;KRAB+ZF nb
# -----------------------------------------------------------------
rule create_global_krabzf_table:
    """
    Creation of global krab and zf result table
    """
    input:
        krab=pathGTDriftGlobalResults + GLOBAL_RESULTS + "table_results/krab_data.csv",
        zf_tabs=expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/tbl/ZF_tabulated",
            accession=ACCESSNB,
        ),
    output:
        pathGTDriftGlobalResults + GLOBAL_RESULTS + "table_results/krabzf_data.csv",
    params:
        path=pathGTDriftData + "genome_assembly"        
    script:
         "python/krabzf.py"    
        
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
            + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/ZF_domains_summary",
            accession=ACCESSNB,
        ),
    output:
        pathGTDriftGlobalResults + GLOBAL_RESULTS + "table_results/zf_count.csv",
    params:
        path=pathGTDriftData + "genome_assembly"      
    script:
         "python/zf_analysis.py"    
                     


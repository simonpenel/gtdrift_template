# RULES FOR PROCESS_STATS_DOMAIN.SMK 
# ==================================

# ------
# concat
# ------
# Concatenate the results of the simple analysis and the paralogs analysis.

rule concat:
    input:
        expand(
            pathGTDriftData
            + "genome_assembly/{{accession}}/analyses/" + GENOME_RESULTS 
            + "summary_hmmsearch_{{accession}}_{domain}.csv",
            domain=DOMAINS_SIMPLE),  
        expand(pathGTDriftData
            + "genome_assembly/{{accession}}/analyses/" + GENOME_RESULTS 
            + "summary_hmmsearch_{{accession}}_{domain}_curated.csv",
            domain=DOMAINS),              
    output:
        pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
            + "whole_summary.csv",
    params:        
        accession=accession_nb,            
    script:
        "../utils/python/concatenate.py"


# ----------------------
# concatenate_assemblies
# ----------------------
# Concatenate the resulats of all assemblies in a global file

rule concatenate_assemblies:
    input:
         # Whole analyse summary
         # ---------------------
         whole_summary=expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
            + "whole_summary.csv",
            accession=ACCESSNB)
    output:
        # Concatenation of assemblies results
        # -----------------------------------
        pathGTDriftGlobalResults
        + GLOBAL_RESULTS + "results.csv"
    script:
        "../utils/python/merge.py"

# -------------
# info_assembly
# -------------
# Get the info associated to an assembly

rule info_assembly:
    input:
        organisms_file=pathGTDriftData + "organisms_data"
    output:
        info_assembly=temp(pathGTDriftGlobalResults + GLOBAL_RESULTS + "info_{accession}.csv")
    params:
        accession=accession_nb
    script:
        "../utils/python/info_assembly.py"       


# ----------------
# list_assemblies
# ----------------
# Concatenate the infos of all assemblies in a  global list

rule list_assemblies:
    input:
        info_assembly=expand(
            pathGTDriftGlobalResults + GLOBAL_RESULTS + "info_{accession}.csv",
            accession=ACCESSNB)
    output:
        # CList of assemblies info
        # ------------------------
        pathGTDriftGlobalResults
        + GLOBAL_RESULTS + "list_assemblies.csv"
    script:
        "../utils/python/list_assemblies.py"
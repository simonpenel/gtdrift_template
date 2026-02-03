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
# concatenate the resulats of all assemblies in a global file

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
        


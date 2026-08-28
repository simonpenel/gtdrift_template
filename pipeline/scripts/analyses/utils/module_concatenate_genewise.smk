rule concat_genewise:
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
            + "whole_summary_genewise.csv",
    params:        
        accession=accession_nb,            
    script:
        "../utils/python/concatenate_genewise.py"
        
#rule concatenate_assemblies:
#    input:
#         # Whole analyse summary
#         # ---------------------
#         whole_summary=expand(
#            pathGTDriftData
#            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
#            + "whole_summary_genewise.csv",
#            accession=ACCESSNB)
#    output:
#        # Concatenation of assemblies results
#        # -----------------------------------
#        pathGTDriftGlobalResults
#        + GLOBAL_RESULTS + "results.csv"
#    script:
#        "../utils/python/merge.py"
        
        
        
rule concatenate_parsed_genwise_results:
    input:
         # Whole analyse summary
         # ---------------------
         whole_summary=expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
            + "parsed_whole_summary_genewise.csv",
            accession=ACCESSNB)
    output:
        # Concatenation of assemblies results
        # -----------------------------------
        temp(pathGTDriftGlobalResults
        + GLOBAL_RESULTS + "parsed_results.csv")
    script:
        "../utils/python/merge.py"
        
        
rule reorder_parsed_genwise_results:
    input:
        # Concatenation of assemblies results
        # -----------------------------------
        pathGTDriftGlobalResults
        + GLOBAL_RESULTS + "parsed_results.csv"
    output:
        # Concatenation of assemblies results
        # -----------------------------------
        pathGTDriftGlobalResults
        + GLOBAL_RESULTS + "ordered_parsed_results.csv"
    script:
        "../utils/python/reorder.py"        
        
rule concatenate_was_parsed_genwise_results:
    input:
         # Whole analyse summary
         # ---------------------
         whole_summary=expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
            + "wad_parsed_whole_summary_genewise.csv",
            accession=ACCESSNB)
    output:
        # Concatenation of assemblies results
        # -----------------------------------
        pathGTDriftGlobalResults
        + GLOBAL_RESULTS + "wad_parsed_results.csv"
    script:
        "../utils/python/merge.py"

# --------------
# info_assembly
# ---------------
# Get the information associated to an assembly

rule info_assembly:
    input:
        organisms_file=pathGTDriftData + "organisms_data"
    output:
        info_assembly=temp(pathGTDriftGlobalResults + GLOBAL_RESULTS + "info_{accession}.csv")
    params:
        accession=accession_nb
    script:
        "../utils/python/info_assembly.py"       

# ----------------------
# list_assemblies
# ----------------------
# Concatenate the infomration on  assemblies in a  list

rule list_assemblies:
    input:
        info_assembly=expand(
            pathGTDriftGlobalResults + GLOBAL_RESULTS + "info_{accession}.csv",
            accession=ACCESSNB)
    output:
        # List of assemblies informations
        # --------------------------------
        pathGTDriftGlobalResults
        + GLOBAL_RESULTS + "list_assemblies.csv"
    script:
        "../utils/python/list_assemblies.py"
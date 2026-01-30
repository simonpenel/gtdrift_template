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
#            + "whole_summary.csv",
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
        pathGTDriftGlobalResults
        + GLOBAL_RESULTS + "parsed_results.csv"
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

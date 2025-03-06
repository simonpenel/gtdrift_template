rule statistics_simple:
    input:
        res=expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
            + "candidates_simple_{{domain}}.txt",accession=ACCESSNB),
        protein=expand(pathGTDriftData + "genome_assembly/{accession}/annotation/protein.faa",accession=ACCESSNB)   
    output:
        pathGTDriftGlobalResults
        + GLOBAL_RESULTS + "simple_statistics_{domain}.csv" 
#    params:    
#        accession=accession_nb,          
    script:
        "../utils/python/statistics_simple.py"

rule simple_statistics_domain:
    input:
        res=pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
            + "candidates_simple_{domain}.txt",
        protein=pathGTDriftData 
            + "genome_assembly/{accession}/annotation/protein.faa"
    output:
           pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
            + "simple_statistics_{domain}.csv"
    params:        
        domain=get_domain, 
        accession=accession_nb,
        type="simple"      
    script:
        "../utils/python/statistics_simple.py"   

rule statistics_domain:
    input:
        res=pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
            + "candidates_{domain}.csv",
        protein=pathGTDriftData 
            + "genome_assembly/{accession}/annotation/protein.faa"
    output:
           pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
            + "statistics_{domain}.csv"
    params:        
        domain=get_domain, 
        accession=accession_nb,
        type="confirmed"     
    script:
        "../utils/python/statistics_simple.py"  

# -------------------------
# statistics_candidates_wad
# -------------------------
# For a genome, give:
# - the number of sequences with all the domains (wad)
# - the number of sequences
# - the taxid and the species 

rule statistics_candidates_wad:
    input:
        res=pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
            + "selected_candidates.csv",
        protein=pathGTDriftData 
            + "genome_assembly/{accession}/annotation/protein.faa",
        organisms_file=pathGTDriftData + "organisms_data"  
    output:
           pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
            + "statistics_wad.csv"
    params:        
        accession=accession_nb,     
    script:
        "../utils/python/statistics_wad.py"  
   
rule statistics_all_domain:
    input:
           domains=expand(pathGTDriftData
            + "genome_assembly/{{accession}}/analyses/" + GENOME_RESULTS 
            + "statistics_{domain}.csv",domain=DOMAINS),
           domains_simple=expand(pathGTDriftData
            + "genome_assembly/{{accession}}/analyses/" + GENOME_RESULTS 
            + "simple_statistics_{domain}.csv",domain=DOMAINS+DOMAINS_SIMPLE),
           organisms_file=pathGTDriftData + "organisms_data"         
    output:
           pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
            + "statistics_all_domains.csv"    
    script:
        "../utils/python/concatenate_domain_stats.py"          

# ----------------------------------
# select_candidates_with_all_domains
# ----------------------------------

rule select_candidates_with_all_domains:
    input:
         candidates_confirmed=expand(pathGTDriftData 
         + "genome_assembly/{{accession}}/analyses/" 
         + GENOME_RESULTS + "candidates_{domain}.csv",
         domain=DOMAINS),
         candidates_simple=expand(pathGTDriftData 
         + "genome_assembly/{{accession}}/analyses/" 
         + GENOME_RESULTS + "candidates_simple_{domain}.txt",
         domain=DOMAINS_SIMPLE),                      
    output:
           pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
            + "selected_candidates.csv"    
    script:
        "../utils/python/select_candidates_wad.py" 
        
rule statistics_candidate_summary:     
    input:
        all_domains_stats=expand(pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
            + "statistics_wad.csv",
            accession=ACCESSNB),
    output:
        pathGTDriftGlobalResults + GLOBAL_RESULTS + "candidate_statistics_summary.csv"
    script:
        "../utils/python/merge_stats.py"                 


rule statistics_all_domain_summary:     
    input:
        all_domains_stats=expand(pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
            + "statistics_all_domains.csv",
            accession=ACCESSNB)
    output:
        pathGTDriftGlobalResults + GLOBAL_RESULTS + "statistics_summary.csv"
    script:
        "../utils/python/merge_stats.py"  
                   
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

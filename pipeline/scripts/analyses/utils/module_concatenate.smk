rule statistics_simple:
    input:
        res=expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
            + "candidates_1_ID_{{domain}}.txt",accession=ACCESSNB),
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
            + "candidates_1_ID_{domain}.txt",
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
    shell:
        """
        cat {input} > {output}.tmp &&
        head -1 {output}.tmp > {output}.tmp2 &&
        grep -v SeqID {output}.tmp >> {output}.tmp2 &&
        cp {output}.tmp2 {output}
        rm {output}.tmp {output}.tmp2 
        """

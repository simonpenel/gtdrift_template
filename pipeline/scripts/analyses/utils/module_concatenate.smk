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
    script:
        "../utils/python/concatenate.py"

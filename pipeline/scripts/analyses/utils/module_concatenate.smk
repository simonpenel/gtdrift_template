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
        + GLOBAL_RESULTS + "results.txt"
    shell:
        """
        cat {input} > {output}.tmp &&
        head -1 {output}.tmp > {output}.tmp2 &&
        grep -v SeqID {output}.tmp >> {output}.tmp2 &&
        cp {output}.tmp2 {output}
        rm {output}.tmp {output}.tmp2 
        """

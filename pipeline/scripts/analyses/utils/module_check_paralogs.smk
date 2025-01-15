import json
import os
import pandas as pd


#rule generate_prdm9_candidates_IDs:
#    input:
#        candidates = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + #"summary_hmmsearch_{accession}.csv",
#    output:
#        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + #"candidates_1_ID.txt"
#    script:
#        "../utils/python/extract_domain_candidates_ID.py"


rule generate_domain_candidates_IDs:
    input:
        domain_per_sequence_tabulated=expand(pathGTDriftData + "genome_assembly/{{accession}}/analyses/" + GENOME_RESULTS + "hmm_search/tbl/{domain}_tabulated",domain=DOMAINS),
        #candidates = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "summary_hmmsearch_{accession}.csv",
    output:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates_1_ID_{domain}.txt"    
    shell:
        """
        cut -f1 {input.domain_per_sequence_tabulated} > {output.candidate_list}
        """
        
        
rule run_seqkit_extract:
    input:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates_1_ID_{domain}.txt",
        multifasta = pathGTDriftData + "genome_assembly/{accession}/annotation/protein.faa"
    output:
        fasta_output = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates_1_{domain}.fasta"
    shell:
        """
        seqkit grep -f {input.candidate_list} {input.multifasta} -o {output.fasta_output}
        """
rule hmm_for_paralogs:
    input:
        alignments_dir = expand(pathGTDriftResource + "ref_align_for_paralogy_check/" + DOMAIN_HMM_DIR +"{domain}",domain=DOMAINS),
    output:
        hmm_db = expand(pathGTDriftResource + "hmm_build/paralogy_check/" + DOMAIN_HMM_DIR +"{domain}.hmm",domain=DOMAINS),
    shell:
        """
        for file in `ls /beegfs/banque/peneldb/gtdrift_template/pipeline/resources/ref_align_for_paralogy_check/domain_example/SET/*fst`; do echo $file; done &

        """       
        
rule hmmscan:
    input:
        hmm_db = pathGTDriftResource + "hmm_build/paralogy_check/" + DOMAIN_HMM_DIR +"profil_{domain}.hmm",
        fasta = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates_1_{domain}.fasta",
        candidate_table = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "summary_hmmsearch_{accession}_{domain}.csv"
    output:
        hmm = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "results_HMM_score_ratio_{domain}.csv",
        table = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "table_score_ratios_HMM_{domain}.csv",
        candidate_table_curated = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "summary_hmmsearch_{accession}_{domain}_curated.csv"
    shell:
        """
        python3 ../utils/python/hmmscan_score_ratio_single.py --hmm_db {input.hmm_db} --output_name {output.hmm} --input_fasta {input.fasta} --output_table {output.table} --csv_file {input.candidate_table}
        """ 

rule curate_prdm9_candidates:
    input:
        candidate_table_curated = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "summary_hmmsearch_{accession}_{domain}_curated.csv",
        assembly_info_file = pathGTDriftData + "organisms_data"
    output:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/"  + GENOME_RESULTS +  "candidates_{domain}.csv"
    params:
        domain_reference=REFERENCE,
    script:
        "../utils/python/select_candidates_bestmatch.py"



rule curate_prdm9_candidates_IDs:
    input:
         candidates = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates_{domain}.csv",
    output:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates_ID_{domain}.txt"
    script:
        "../utils/python/extract_domain_candidates_ID.py"

rule run_seqkit_extraction_curated:
    input:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates_ID_{domain}.txt",
        multifasta = pathGTDriftData + "genome_assembly/{accession}/annotation/protein.faa"
    output:
        fasta_output = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates_{domain}.fasta"
    shell:
        """
        seqkit grep -f {input.candidate_list} {input.multifasta} -o {output.fasta_output}
        """







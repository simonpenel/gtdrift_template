import json
import os
import pandas as pd


rule generate_prdm9_candidates_IDs:
    input:
        candidates = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "summary_hmmsearch_{accession}.csv",
    output:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates_1_ID.txt"
    script:
        "../utils/python/extract_domain_candidates_ID.py"

rule run_seqkit_extract:
    input:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates_1_ID.txt",
        multifasta = pathGTDriftData + "genome_assembly/{accession}/annotation/protein.faa"
    output:
        fasta_output = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates_1.fasta"
    shell:
        """
        seqkit grep -f {input.candidate_list} {input.multifasta} -o {output.fasta_output}
        """

rule hmmscan:
    input:
        hmm_db = pathGTDriftResource + "ref_align/" + REFERENCE_DOMAIN_HMM,
        fasta = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates_1.fasta",
        candidate_table = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "summary_hmmsearch_{accession}.csv"
    output:
        hmm = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "results_HMM_score_ratio.csv",
        table = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "table_score_ratios_HMM.csv",
        candidate_table_curated = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "summary_hmmsearch_{accession}_curated.csv"
    shell:
        """
        python3 ../utils/python/hmmscan_score_ratio.py --hmm_db {input.hmm_db} --output_name {output.hmm} --input_fasta {input.fasta} --output_table {output.table} --csv_file {input.candidate_table}
        """ 

rule curate_prdm9_candidates:
    input:
        candidate_table_curated = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "summary_hmmsearch_{accession}_curated.csv",
        assembly_info_file = pathGTDriftData + "organisms_data"
    output:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/"  + GENOME_RESULTS +  "candidates.csv"
    params:
        domain_reference=REFERENCE,
    script:
        "../utils/python/select_candidates_bestmatch.py"
#    shell:
#       """
#       python3 ../utils/python/select_candidates_bestmatch.py {input.candidate_table_curated} {input.assembly_info_file} {output}
#       """



rule curate_prdm9_candidates_IDs:
    input:
         candidates = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates.csv",
    output:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates_ID.txt"
    script:
        "../utils/python/extract_domain_candidates_ID.py"

rule run_seqkit_extraction_curated:
    input:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates_ID.txt",
        multifasta = pathGTDriftData + "genome_assembly/{accession}/annotation/protein.faa"
    output:
        fasta_output = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates.fasta"
    shell:
        """
        seqkit grep -f {input.candidate_list} {input.multifasta} -o {output.fasta_output}
        """







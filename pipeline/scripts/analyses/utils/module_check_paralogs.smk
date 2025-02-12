import json
import os
import pandas as pd

# ---------------------------------------------------------------------------------------
# generate_domain_candidates_IDs
# Extract the list of candidates from the search of a hmm profile domain in the proteome.
# ---------------------------------------------------------------------------------------
rule generate_domain_candidates_IDs:
    input:
        domain_per_sequence_tabulated=pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/tbl/{domain}_tabulated"
    output:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates_1_ID_{domain}.txt"    
    shell:
        """
        echo process {input.domain_per_sequence_tabulated} &&
        cut -f1 {input.domain_per_sequence_tabulated} > {output.candidate_list}
        """
        
# --------------------------------------------------------------------
# run_seqkit_extract
# Extract the fasta sequence of the candidates from its sequence name.
# --------------------------------------------------------------------      
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
 
# -------------------------------------------------------------
# hmmscan
# Scan the candidates into the hmm database and calculate ratio
# of the best match score according to the next match.
# -------------------------------------------------------------      
rule hmmscan:
    input:
        hmm_db = pathGTDriftResource + RESOURCES_DIR_NAME + "hmm_databases/database_{domain}.hmm",
        hmm_database_index1=pathGTDriftResource + RESOURCES_DIR_NAME + "hmm_databases/database_{domain}.hmm.h3f", 
        hmm_database_index2=pathGTDriftResource + RESOURCES_DIR_NAME + "hmm_databases/database_{domain}.hmm.h3i",
        hmm_database_index3=pathGTDriftResource + RESOURCES_DIR_NAME + "hmm_databases/database_{domain}.hmm.h3m",
        hmm_database_index4=pathGTDriftResource + RESOURCES_DIR_NAME + "hmm_databases/database_{domain}.hmm.h3p",  
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

# -----------------------------------------------------------
# curate_prdm9_candidates
# Select the candidate if the best hmm match is the reference.
# ------------------------------------------------------------
rule curate_prdm9_candidates:
    input:
        candidate_table_curated = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "summary_hmmsearch_{accession}_{domain}_curated.csv",
        assembly_info_file = pathGTDriftData + "organisms_data"
    output:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/"  + GENOME_RESULTS +  "candidates_{domain}.csv"
    params:
        domain_reference=get_reference,
        domain=get_domain,
    script:
        "../utils/python/select_candidates_bestmatch.py"

# --------------------------------
# curate_prdm9_candidates_IDs
# Get the selected candidates ids.
# --------------------------------
rule curate_prdm9_candidates_IDs:
    input:
         candidates = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates_{domain}.csv",
    output:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates_ID_{domain}.txt"
    script:
        "../utils/python/extract_domain_candidates_ID.py"


# -----------------------------------------------------------------------------
# run_seqkit_extraction_curated
# Extract the fasta sequence of the selected candidates from its sequence name.
# ------------------------------------------------------------------------------  
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


# ---------------------------------------------------------------------------------------
# generate_domain_candidates_wad_IDs
# Extract the list of candidates with all domains from the csv from the search of a hmm profile domain in the proteome.
# ---------------------------------------------------------------------------------------
rule generate_candidates_wad_IDs:
    input:
        candidates = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates_with_all_domains.csv"
    output:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates_with_all_domains.txt"   
    script:
        "../utils/python/extract_domain_candidates_ID.py"
        
        
        
# -----------------------------------------------------------------------------
# run_seqkit_extraction_candidates_wad
# Extract the fasta sequence of the candidates from its sequence name.
# ------------------------------------------------------------------------------  
rule run_seqkit_extraction_wad_fasta:
    input:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates_with_all_domains.txt",
        multifasta = pathGTDriftData + "genome_assembly/{accession}/annotation/protein.faa"
    output:
        fasta_output = pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "candidates.fasta"
    shell:
        """
        seqkit grep -f {input.candidate_list} {input.multifasta} -o {output.fasta_output}
        """










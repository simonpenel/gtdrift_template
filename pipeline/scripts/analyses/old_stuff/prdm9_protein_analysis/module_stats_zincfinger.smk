import json
import os
import pandas as pd

configfile: "config.json"
ACCESSNB = config["assembly_list"]

# Function to load JSON files
def load_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

# Assign environment variables
globals().update(load_json("../environment_path.json"))


rule generate_prdm9_candidates:
    input:
        prdm9_prot_summary = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/summary_hmmsearch_prdm9_with_paralog_check_{accession}.csv",
        assembly_info_file = pathGTDriftData + "organisms_data"
    output:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/prdm9_candidates_1.csv"
    shell:
        """
        python3 python/generate_PRDM9_candidates.py {input.prdm9_prot_summary} {input.assembly_info_file} {output}
        """

rule generate_prdm9_candidates_IDs:
    input:
        prdm9_prot_dir = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/prdm9_candidates_1.csv"
    output:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/prdm9_candidates_1_ID.txt"
    shell:
        """
        python3 python/extract_PRDM9_candidates_ID.py {input} {output}
        """

rule run_seqkit_extract:
    input:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/prdm9_candidates_1_ID.txt",
        multifasta = pathGTDriftData + "genome_assembly/{accession}/annotation/protein.faa"
    output:
        fasta_output = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/candidates_prdm9_1.fasta"
    shell:
        """
        seqkit grep -f {input.candidate_list} {input.multifasta} -o {output.fasta_output}
        """

rule hmmscan:
    input:
        hmm_db = pathGTDriftResource + "ref_align/SET_PRDMs_Metazoa_Reference_alignment/concatenated_SET_PRDMs_alignments.hmm",
        fasta = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/candidates_prdm9_1.fasta",
        candidate_table = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/summary_hmmsearch_prdm9_{accession}.csv"
    output:
        hmm = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/results_HMM_score_ratio.csv",
        table = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/table_score_ratios_HMM.csv",
        candidate_table_curated = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/summary_hmmsearch_prdm9_{accession}_curated.csv"
    shell:
        """
        python3 python/hmmscan_score_ratio.py --hmm_db {input.hmm_db} --output_name {output.hmm} --input_fasta {input.fasta} --output_table {output.table} --csv_file {input.candidate_table}
        """ 

rule curate_prdm9_candidates:
    input:
        candidate_table_curated = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/summary_hmmsearch_prdm9_{accession}_curated.csv",
        assembly_info_file = pathGTDriftData + "organisms_data"
    output:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/prdm9_candidates.csv"
    shell:
        """
        python3 python/generate_PRDM9_candidates.py {input.candidate_table_curated} {input.assembly_info_file} {output}
        """

rule curate_prdm9_candidates_IDs:
    input:
        prdm9_prot_dir = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/prdm9_candidates.csv"
    output:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/prdm9_candidates_ID.txt"
    shell:
        """
        python3 python/extract_PRDM9_candidates_ID.py {input} {output}
        """

rule run_seqkit_extraction_curated:
    input:
        candidate_list = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/prdm9_candidates_ID.txt",
        multifasta = pathGTDriftData + "genome_assembly/{accession}/annotation/protein.faa"
    output:
        fasta_output = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/candidates_prdm9.fasta"
    shell:
        """
        seqkit grep -f {input.candidate_list} {input.multifasta} -o {output.fasta_output}
        """

rule merge_prdm9_candidates:
    """
    Merge all individual prdm9_candidates.csv files into a single global file.
    """
    input:
        expand(pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/prdm9_candidates.csv", accession=ACCESSNB)
    output:
        merged_candidates = pathGTDriftGlobalResults + GLOBAL_RESULTS + "table_results/global_prdm9_candidates.csv"
    run:
        # Combine all CSV files into one DataFrame
        dfs = [pd.read_csv(file, sep=';') for file in input]
        combined_df = pd.concat(dfs, ignore_index=True)
        combined_df.to_csv(output.merged_candidates, sep=';', index=False)

rule zincfinger_analysis:
    """
    Run the zinc finger analysis on each protein sequence using R.
    """
    input:
        protein_seq = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/candidates_prdm9.fasta"
    output:
        zincfinger_out = pathGTDriftGlobalResults + GLOBAL_RESULTS + "zinc_finger/ZFD_{accession}.csv"
    
    run:
        # Verify if the input file exists
        if not os.path.exists(input.protein_seq):
            raise FileNotFoundError(f"Input file does not exist: {input.protein_seq}")
        
        # Run the R script
        command = f"Rscript --vanilla ./zincfinger_analysis.R {input.protein_seq} {output.zincfinger_out}"
        shell(command)

rule combine_zinc_finger:
    """
    Combine all ZFD_{accession}.csv files into one zinc_finger.csv.
    """
    input:
        zfd_files = expand(pathGTDriftGlobalResults + GLOBAL_RESULTS + "zinc_finger/ZFD_{accession}.csv", accession=ACCESSNB)
    output:
        zinc_finger_combined = pathGTDriftGlobalResults + GLOBAL_RESULTS + "table_results/zinc_finger.csv"
    run:
        # Check if there are input files before trying to combine
        if not input.zfd_files:
            raise FileNotFoundError("No input files found for combining.")
        
        # Combine all ZFD_{accession}.csv files into one zinc_finger.csv
        with open(output.zinc_finger_combined, 'w') as outfile:
            for fname in input.zfd_files:
                with open(fname) as infile:
                    outfile.write(infile.read())
                    outfile.write("\n")  # Add a newline between files

rule general_table:
    """
    Generate a general table of PRDM9 candidates
    """
    input:
        merged_candidates = pathGTDriftGlobalResults + GLOBAL_RESULTS + "table_results/global_prdm9_candidates.csv"
    output:
        general_table = pathGTDriftGlobalResults + GLOBAL_RESULTS + "table_results/table_prdm9.csv"
    shell:
        """
        python3 python/general_table_prdm9.py -i {input} -o {output}
        """


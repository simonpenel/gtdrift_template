import pandas as pd
import argparse

'''
This script takes as an input the output .csv file from the script "generate_PRDM9_candidates.py" and generates
a .txt file with the IDs of all PRDM9 protein candidates
'''

def extract_seq_ids(input_file, output_file):
    try:
        # Attempt to read the CSV file
        df = pd.read_csv(input_file, sep=';')
    except pd.errors.EmptyDataError:
        # If the file is empty, create an empty output file
        with open(output_file, 'w') as f:
            pass
        return

    # Check if the DataFrame is empty or if it lacks the 'SeqID' column
    if df.empty or 'SeqID' not in df.columns:
        # If the file is empty or missing the required column, create an empty file
        with open(output_file, 'w') as f:
            pass
    else:
        # Extract the 'SeqID' column
        seq_ids = df['SeqID']

        # Write each `SeqID` on a new line in the output file
        with open(output_file, 'w') as f:
            for seq_id in seq_ids:
                f.write(f"{seq_id}\n")

if __name__ == "__main__":

    # Call the function to extract SeqIDs
    extract_seq_ids(snakemake.input.candidates, snakemake.output.candidate_list)

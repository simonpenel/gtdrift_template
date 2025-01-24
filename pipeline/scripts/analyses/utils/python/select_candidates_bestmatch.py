import pandas as pd
import argparse
import os

'''
This script takes the "summary_table_(Genome_Accesion).csv" file as input, filters only

'''

def filter_bestmatch(input_file, output_file, reference, domain):
    try:
        # Read the original CSV file using ';' as the separator
        df = pd.read_csv(input_file, sep=';', header=0)
    except pd.errors.EmptyDataError:
        # If the input file is empty, generate the default output row
        create_default_output(output_file, domain)
        return

    # Remove the 'Unnamed' columns if they exist
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

    # Check if the DataFrame is empty after reading
    if df.empty:
        # If only the headers are present, generate the default output row
        #create_default(output_file)
        create_default_output(output_file, domain)
    else:
        # Filter the rows where the 'Best Match' column equals 'PRDM9'
        df_filtered = df[df['Best Match'] == reference ]
        
        # If the filtered DataFrame has no rows, generate a row with default values
        if df_filtered.empty:
            #create_default(output_file)
            create_default_output(output_file, domain)
            return
        

        # Save the resulting DataFrame to a new CSV file with semicolon delimiter
        df_filtered.to_csv(output_file, sep=';', index=False)


def create_default_output(output_file, domain):
    # Create the output file when no valid rows are found in the input file or if the file is empty
    # Generate a row with empty values for the columns
    new_row = {
        'SeqID': "",
        domain + ' Query': "",
        domain + ' E-value': "",
        domain + ' Score': "",
        'Nb '+domain + ' domains': "",
        domain + ' domain start': "",
        domain + ' domain end': "",
        'Taxid': "",
        'Best Match': "",
        'Bit Score': "",
        'Score ratio': ""
    }
    # Convert the dictionary to a DataFrame and add it to the output DataFrame
    new_df = pd.DataFrame([new_row])
    # Save the output DataFrame with the generated line
    new_df.to_csv(output_file, sep=';', index=False)
    
if __name__ == "__main__":
    domain = snakemake.params.domain
    domain_reference = snakemake.params.domain_reference
    candidates_file = snakemake.input.candidate_table_curated
    selected_candidates_list_file = snakemake.output.candidate_list
    # Call the function to perform the filtering
    filter_bestmatch(candidates_file, selected_candidates_list_file, domain_reference, domain)

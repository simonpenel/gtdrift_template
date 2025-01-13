import pandas as pd
import argparse
import os

'''
This script takes the "summary_table_(Genome_Accesion).csv" file as input, filters only

'''

def filter_bestmatch(input_file, output_file,reference):
    # Get the absolute path of the input file and move up three levels to get the genome assembly directory
    #input_dir = os.path.dirname(os.path.abspath(input_file))  # Get the directory of the input file
    #genome_assembly = os.path.basename(os.path.dirname(os.path.dirname(input_dir)))  # Move up three levels

    try:
        # Read the original CSV file using ';' as the separator
        df = pd.read_csv(input_file, sep=';', header=0)
    except pd.errors.EmptyDataError:
        # If the input file is empty, generate the default output row
        create_default(output_file)
        #create_default_output(output_file, genome_assembly, taxid_species_file)
        return

    # Remove the 'Unnamed' columns if they exist
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

    # Check if the DataFrame is empty after reading
    if df.empty:
        # If only the headers are present, generate the default output row
        create_default(output_file)
    else:
        # Filter the rows where the 'Best Match' column equals 'PRDM9'
        df_filtered = df[df['Best Match'] == reference ]
        
        # If the filtered DataFrame has no rows, generate a row with default values
        if df_filtered.empty:
            create_default(output_file)
            return
        

        # Save the resulting DataFrame to a new CSV file with semicolon delimiter
        df_filtered.to_csv(output_file, sep=';', index=False)

def create_default(output_file):
    print("\n\n\nNO DATA\n\n\n")

def create_default_output(output_file, genome_assembly, taxid_species_file):
    # Create the output file when no valid rows are found in the input file or if the file is empty
    try:
        taxid_species_df = pd.read_csv(taxid_species_file, sep='\t', header=0)
    except pd.errors.EmptyDataError:
        # If the new file is empty, create an empty output file
        pd.DataFrame(columns=['species_name', 'taxid', 'Assembly']).to_csv(output_file, sep=';', index=False)
        return

    # Look for the value of 'species_name' and 'taxid' using the 'Assembly'
    species_names = []
    taxids = []
    for assembly in [genome_assembly]:
        # Look in taxid_species_df where 'Assembly Accession' matches and extract 'Species Name' and 'Taxid'
        match = taxid_species_df[taxid_species_df['Assembly Accession'] == assembly]
        if not match.empty:
            species_names.append(match['Species Name'].values[0])
            taxids.append(match['Taxid'].values[0])
        else:
            species_names.append(None)  # If not found, leave None
            taxids.append(None)         # If not found, leave None

    # Generate a row with empty values for the columns
    new_row = {
        'species_name': species_names[0],
        'taxid': taxids[0],
        'Assembly': genome_assembly,
        'SeqID': "",
        'SET Query': "",
        'SET E-value': "",
        'SET Score': "",
        'Nb SET domains': 0,
        'SET domain start': "",
        'SET domain end': "",
        'KRAB Query': "",
        'KRAB E-value': "",
        'KRAB Score': "",
        'Nb KRAB domains': 0,
        'KRAB domain start': "",
        'KRAB domain end': "",
        'SSXRD Query': "",
        'SSXRD E-value': "",
        'SSXRD Score': "",
        'Nb SSXRD domains': 0,
        'SSXRD domain start': "",
        'SSXRD domain end': "",
        'ZF Query': "",
        'ZF E-value': "",
        'ZF Score': "",
        'Nb ZF domains': 0,
        'ZF domain start': "",
        'ZF domain end': "",
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
    # Create the argument parser
    #parser = argparse.ArgumentParser(description="Filter sequences where the best match is PRDM9 and add species information.")

    # Add the arguments for the input file, taxid/species file, and output file
    #parser.add_argument('input_file', type=str, help="Path to the input CSV file")
    #parser.add_argument('taxid_species_file', type=str, help="Path to the TSV file with taxid and species_name information")
    #parser.add_argument('output_file', type=str, help="Path to the output file (filtered result)")

    # Parse the arguments
    #args = parser.parse_args()
    domain_reference = snakemake.params.domain_reference
    candidates_file = snakemake.input.candidate_table_curated
    #organisms_file = snakemake.input.organisms_file
    selected_candidates_list_file = snakemake.output.candidate_list
    # Call the function to perform the filtering
    filter_bestmatch(candidates_file, selected_candidates_list_file,domain_reference)

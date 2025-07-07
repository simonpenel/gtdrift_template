import sys
import csv
import os
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

# This python script identifies in the proteic sequence/s of a FASTA the three positions in which it should be the three 
# catalytic tyrosines of the SET domain of the PRDM9 gene and checks if there are tyrosines there (1) or not (0)

# Main function
def run_blastp(query_seq, subject_file, output_csv):
    # Write the query sequence to a temporary file
    query_file = "query.fasta"
    with open(query_file, "w") as f:
        f.write(f">Query\n{query_seq}\n")
    
    # Run blastp
    blastp_cline = NcbiblastpCommandline(query=query_file, subject=subject_file, outfmt=5, out="blast_result.xml")
    stdout, stderr = blastp_cline()

    # Parse the output XML file
    with open("blast_result.xml") as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        
        # Prepare a dictionary to store the best hits for each subject
        best_hits = {}

        # Process each BLAST hit
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    subject_id = alignment.hit_id
                    print("debug sunject_id = "+subject_id)
                    # Check if this subject_id has been seen before
                    if subject_id not in best_hits:
                        best_hits[subject_id] = hsp  # Store HSP as the best hit for this subject
                    else:
                        # Compare the current HSP score with the stored best hit score
                        print("debug "+str(hsp.score) + "> " + str(best_hits[subject_id].score))
                        if hsp.score > best_hits[subject_id].score:
                            best_hits[subject_id] = hsp  # Update the best hit

        # Prepare the output CSV file
        with open(output_csv, "w", newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            # Write the header with new column names
            csvwriter.writerow(["Subject_ID", "Y276", "Y341", "Y357"])

            # Write the best hits to the CSV
            for subject_id, hsp in best_hits.items():
                # Get aligned portions of query and subject
                aligned_query = hsp.query[0:hsp.align_length]  # Portion of the query aligned
                aligned_subject = hsp.sbjct[0:hsp.align_length]  # Portion of the subject aligned

                # Initialize amino acid indicators for the subject
                y276 = y341 = y357 = 0

                # Variables to track the current positions in the aligned sequences
                query_position_in_alignment = hsp.query_start
                subject_position_in_alignment = hsp.sbjct_start

                # Variables to track the real positions in the sequences (ignoring gaps from indels introduced during the blastp)
                query_real_pos = query_position_in_alignment
                subject_real_pos = subject_position_in_alignment

                # Iterate over the aligned query and subject sequences
                for idx, (query_char, subject_char) in enumerate(zip(aligned_query, aligned_subject)):
                    if query_char != '-':  # Ignore gaps in query
                        # Check if we are at the positions of interest in the query
                        if query_real_pos == 72:
                            y276 = 1 if subject_char == 'Y' else 0  # Corresponding to Y276 
                        elif query_real_pos == 137:
                            y341 = 1 if subject_char == 'Y' else 0  # Corresponding to Y341
                        elif query_real_pos == 153:
                            y357 = 1 if subject_char == 'Y' else 0  # Corresponding to Y357

                        query_real_pos += 1  # Move to the next real position in the query sequence

                    if subject_char != '-':  # Ignore gaps in subject
                        subject_real_pos += 1  # Move to the next real position in the subject sequence

                # Write the data to the CSV with updated column names
                csvwriter.writerow([subject_id, y276, y341, y357])

def is_fasta_empty(subject_file):
    """ Check if the FASTA file is empty or does not exist. """
    return not os.path.exists(subject_file) or os.path.getsize(subject_file) == 0


# Check the arguments passed from the command line
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python blastp_script.py <subject_file> <output_csv>")
        sys.exit(1)

    # Reference (query sequence), in this case the SET domain of the PRDM9 human gene
    query_seq = "CEMCQNFFIDSCAAHGPPTFVKDSAVDKGHPNRSALSLPPGLRIGPSGIPQAGLGVWNEASDLPLGLHFGPYEGRITEDEEAANNGYSWLITKGRNCYEYVDGKDKSWANWMRYVNCARDDEEQNLVAFQYHRQIFYRTCRVIRPGCELLVWYGDEYGQELGIKWGSKWKKELMAGR"
    
    # Subject sequences file passed as an argument
    subject_file = sys.argv[1]
    
    # Output CSV file
    output_csv = sys.argv[2]
    
# Check if the FASTA file is empty
    if is_fasta_empty(subject_file):
        # Generate output CSV with the message
        with open(output_csv, "w", newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(["0 candidate PRDM9 proteins"])
    else:
        # Run BLASTP and save the results in the CSV file
        run_blastp(query_seq, subject_file, output_csv)


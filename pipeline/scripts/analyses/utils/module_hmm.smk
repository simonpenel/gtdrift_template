## Main snakemake for domain analysis on protein data
## Date : Decembre 2024
## Authors :


# Import python modules
# ---------------------
import os
import json

#directories, files = glob_wildcards(pathGTDriftResource + ANALYSE_DIR + "reference_alignments/{dir}/{file}.fst") 

# Rules
# -----

# --------------------------------------------------------------
# all : inputs define the files generated at the end of process. 
# --------------------------------------------------------------
##rule all_bidon:
#    """
#    Get the candidates in fasta format
#    """
#    input:       
#        # all hmm      
#        hmm_profiles=expand(pathGTDriftResource + ANALYSE_DIR + "hmm_profiles/{dir}/{file}.hmm", #dir=directories, file=files),
#        hmm_database=expand(pathGTDriftResource + ANALYSE_DIR + "hmm_databases/database_{dir}.hmm", #dir=directories),
#        hmm_database_index1=expand(pathGTDriftResource + ANALYSE_DIR + "hmm_databases/database_{dir}.hmm.h3f", #dir=directories),   
#        hmm_database_index2=expand(pathGTDriftResource + ANALYSE_DIR + "hmm_databases/database_{dir}.hmm.h3i", #dir=directories),        
#        hmm_database_index3=expand(pathGTDriftResource + ANALYSE_DIR + "hmm_databases/database_{dir}.hmm.h3m", #dir=directories),
#        hmm_database_index4=expand(pathGTDriftResource + ANALYSE_DIR + "hmm_databases/database_{dir}.hmm.h3p", #dir=directories),           

  
rule calc_hmm_profiles:
    input:
        pathGTDriftResource + ANALYSE_DIR + "reference_alignments/{dir}/{file}.fst",
    output:
        pathGTDriftResource + ANALYSE_DIR + "hmm_profiles/{dir}/{file}.hmm",
    shell:
        "{RUNCMD} hmmbuild {output} {input}"        

rule calc_hmm_database:
    input:
        pathGTDriftResource + ANALYSE_DIR + "hmm_profiles/{dir}",
    output:
        pathGTDriftResource + ANALYSE_DIR + "hmm_databases/database_{dir}.hmm",
    shell:
        "cat {input}/*.hmm > {output}"          

rule index_hmm_database:
    input:
        pathGTDriftResource + ANALYSE_DIR + "hmm_databases/database_{dir}.hmm",
    output:
        pathGTDriftResource + ANALYSE_DIR + "hmm_databases/database_{dir}.hmm.h3f",
        pathGTDriftResource + ANALYSE_DIR + "hmm_databases/database_{dir}.hmm.h3i",
        pathGTDriftResource + ANALYSE_DIR + "hmm_databases/database_{dir}.hmm.h3m",
        pathGTDriftResource + ANALYSE_DIR + "hmm_databases/database_{dir}.hmm.h3p",
    shell:
        "{RUNCMD} hmmpress  {input}"   


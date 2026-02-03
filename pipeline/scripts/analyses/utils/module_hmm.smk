## Module for hmm profile calculation
## Date : Decembre 2024

# Import python modules
# ---------------------
import os
import json

# Rules
# -----

# ------------------------------------------
# calc_hmm_profiles
# Calculate a hmm profile from an alignment.   
# ------------------------------------------  
rule calc_hmm_profiles:
    input:
        pathGTDriftResource + RESOURCES_DIR_NAME + "reference_alignments/{dir}/{file}.fst",
    output:
        pathGTDriftResource + RESOURCES_DIR_NAME + "hmm_profiles/{dir}/{file}.hmm",
    shell:
        "{RUNCMD} hmmbuild {output} {input}"        

# ------------------------------------------------------------
# calc_hmm_database
# Concatenate the hmm profiles of a domain into a hmm database.
# -------------------------------------------------------------
rule calc_hmm_database:
    input:
        dir=pathGTDriftResource + RESOURCES_DIR_NAME + "hmm_profiles/{dir}",
        # This is not an input for the shell, but it should be present to force the calculation
        # of all the hmm profiles if necessary (and not only the reference hmm profile)
        hmm_profiles=expand(pathGTDriftResource + RESOURCES_DIR_NAME + "hmm_profiles/{dir}/{file}.hmm",
               zip, dir=directories, file=files),
    output:
        database=pathGTDriftResource + RESOURCES_DIR_NAME + "hmm_databases/database_{dir}.hmm",
    shell:
        "cat {input.dir}/*.hmm > {output.database}"          

# -------------------------------
# index_hmm_database
# indexation of the hmm database.
# ------------------------------- 
rule index_hmm_database:
    input:
        pathGTDriftResource + RESOURCES_DIR_NAME + "hmm_databases/database_{dir}.hmm",
    output:
        pathGTDriftResource + RESOURCES_DIR_NAME + "hmm_databases/database_{dir}.hmm.h3f",
        pathGTDriftResource + RESOURCES_DIR_NAME + "hmm_databases/database_{dir}.hmm.h3i",
        pathGTDriftResource + RESOURCES_DIR_NAME + "hmm_databases/database_{dir}.hmm.h3m",
        pathGTDriftResource + RESOURCES_DIR_NAME + "hmm_databases/database_{dir}.hmm.h3p",
    shell:
        "{RUNCMD} hmmpress  {input}"   


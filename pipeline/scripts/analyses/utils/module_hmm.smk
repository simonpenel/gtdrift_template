## Main snakemake for domain analysis on protein data
## Date : Decembre 2024
## Authors :


# Import python modules
# ---------------------
import os
import json

# Rules
# -----

           
rule calc_hmm_profiles:
    input:
        pathGTDriftResource + ANALYSE_DIR + "reference_alignments/{dir}/{file}.fst",
    output:
        pathGTDriftResource + ANALYSE_DIR + "hmm_profiles/{dir}/{file}.hmm",
    shell:
        "{RUNCMD} hmmbuild {output} {input}"        

rule calc_hmm_database:
    input:
        dir=pathGTDriftResource + ANALYSE_DIR + "hmm_profiles/{dir}",
        # This is not an input for the shell, but it should be present to force the calculation
        # of all the hmm profiles if necessary (and not only the reference hmm profile)
        hmm_profiles=expand(pathGTDriftResource + ANALYSE_DIR + "hmm_profiles/{dir}/{file}.hmm",
               zip, dir=directories, file=files),
    output:
        database=pathGTDriftResource + ANALYSE_DIR + "hmm_databases/database_{dir}.hmm",
    shell:
        "cat {input.dir}/*.hmm > {output.database}"          

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


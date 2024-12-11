## Main snakemake for prdm9 analysis on protein data
## Date : Decembre 2024
## Authors :


# Import python modules
# ---------------------
import os
import json


# Function to load JSON files
# ---------------------------
def load_json(file_path):
    with open(file_path, "r") as file:
        return json.load(file)


# Assign environment variables
# ----------------------------
globals().update(load_json("../environment_path.json"))


# Configuration
# -------------
configfile: "config.json"


ACCESSNB = config["assembly_list"]
DOMAIN = ["KRAB", "SET", "SSXRD", "ZF"]

# Rules
# -----

# --------------------------------------------------------------
# all : inputs define the files generated at the end of process. 
# --------------------------------------------------------------
rule all:
    """
    Get the prdm9 stats on protein data
    """
    input:
        # prdm9 protein statistics for each assembly.
        # warning : summary_hmmsearch_prdm9_{accession}.csv should
        # no be confused with summary_hmmsearch_prdm9_with_paralog_check_{accession}.csv
        # -----------------------------------------------------
        stats_prdm9=expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/prdm9_prot/summary_hmmsearch_prdm9_{accession}.csv",
            accession=ACCESSNB,
        ),
        # global statistics on KRAB domains
        # ---------------------------------
        krab=pathGTDriftGlobalResults + "analyses_summaries/table_results/krab_data.csv",
        # global statistics on KRAB and ZF domains
        # ----------------------------------------
        ##krabzf=pathGTDriftGlobalResults
        + "analyses_summaries/table_results/krabzf_data.csv",
        # global statistics on ZF domains count
        # -------------------------------------
        ##zf=pathGTDriftGlobalResults + "analyses_summaries/table_results/zf_count.csv",
        # global statistics on prdm9
        # --------------------------
        table=pathGTDriftGlobalResults
        + "analyses_summaries/table_results/table_prdm9.csv",
        # summary of prdm9 candidates
        # ---------------------------
        PRDM9_candidates=pathGTDriftGlobalResults
        + "analyses_summaries/table_results/global_prdm9_candidates.csv",
        # statistics on zinc finger diversity
        # -----------------------------------
        ##zincfinger=pathGTDriftGlobalResults
        + "analyses_summaries/table_results/zinc_finger.csv",
        # statistics on SET tyrosines
        # ---------------------------
        ##SET_tyrosines=pathGTDriftGlobalResults
        + "analyses_summaries/table_results/SET_tyrosines.csv",


# Modules snakemake
# -----------------

# -----------------------------------------------------------
# module_stats_prdm9.smk:
# general statistics on prdm9 search in protein data.
# outputs : stats_prdm9
# -----------------------------------------------------------
include: "module_stats_prdm9.smk"
# -----------------------------------------------------------
# module_stats_zinfinger.smk:
# statistics on zinc fingers in protein data, summary of data
# outputs : PRDM9_candidates, zincfinger, table
# -----------------------------------------------------------
include: "module_stats_zincfinger.smk"
# -----------------------------------------------------------
# module_SET_tyrosines.smk:
# statistics on SET tyrosines
# outputs : SET_tyrosines
# -----------------------------------------------------------
include: "module_SET_tyrosines.smk"
# -----------------------------------------------------------
# module_concatenate_results.smk:
# concatenate
# outputs : krab,krabzf,zf
# -----------------------------------------------------------
include: "module_concatenate_results.smk"

import os

import json

# Function to load JSON files
def load_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

# Assign environment variables
globals().update(load_json("../environment_path.json"))


configfile: "config.json"
ACCESSNB  = config["assembly_list"]
DOMAIN = ['KRAB', 'SET', 'SSXRD', 'ZF']

rule all:
    """
    Get the prdm9 stats
    """
    input:
        stats_prdm9 = expand(pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/summary_table_prdm9_{accession}.csv", accession=ACCESSNB),
        krab = pathGTDriftGlobalResults +"analyses_summaries/table_results/krab_data.csv",
        krabzf = pathGTDriftGlobalResults + "analyses_summaries/table_results/krabzf_data.csv",
        zf = pathGTDriftGlobalResults + "analyses_summaries/table_results/zf_count.csv",
        table = pathGTDriftGlobalResults + "analyses_summaries/table_results/table_prdm9.csv",
        PRDM9_candidates = pathGTDriftGlobalResults + "analyses_summaries/table_results/global_prdm9_candidates.csv",
	zincfinger = pathGTDriftGlobalResults + "analyses_summaries/table_results/zinc_finger.csv",
	SET_tyrosines = pathGTDriftGlobalResults + "analyses_summaries/table_results/SET_tyrosines.csv"

include: "module_stats_prdm9.smk"
include: "module_stats_zincfinger.smk"
include: "module_SET_tyrosines.smk"

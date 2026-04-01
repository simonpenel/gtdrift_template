import pandas as pd
import os
import json


#configfile:"config.json"

# Function to load JSON files
def load_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

# Assign environment variables
globals().update(load_json("../environment_path.json"))

localrules: download_NCBI_genome, download_NCBI_protein, download_NCBI_cds, download_NCBI_annotation

if "assembly_list" in config.keys():
    assembly_list = config["assembly_list"]
else:
    assembly_list = None

storagetype  = config["storagetype"]

rule collect_everything:
     input:
         #expand( pathGTDriftData + "genome_assembly_ebi/{genome_assembly}/genome_seq/genomic.fna.path",genome_assembly=assembly_list),
         expand( pathGTDriftData + "genome_assembly_ebi/{genome_assembly}/annotation/protein.faa",genome_assembly=assembly_list),
         #expand( pathGTDriftData + "genome_assembly_ebi/{genome_assembly}/annotation/cds_from_genomic.fna",genome_assembly=assembly_list),
         #expand( pathGTDriftData + "genome_assembly_ebi/{genome_assembly}/annotation/genomic.gff",genome_assembly=assembly_list)


# Ne se lance que sur un noeud (-j = 1) sinon bug
rule download_NCBI_protein:
    # Telecharge le fichier de sequence proteique de EBI
    params:
        assembly = "{genome_assembly}"
    input:
        data_organisms = "organisms_data_ebi"
    output:
        prot_path =  pathGTDriftData + "genome_assembly_ebi/{genome_assembly}/annotation/protein.faa"
    shell:
        "python3 {pathGTDriftScripts}analyses/collecting_genome_annotation_ebi/download_proteins_ebi.py -d {input.data_organisms} --assembly {params.assembly} --output {output.prot_path}"


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
         expand( pathGTDriftData + "genome_assembly/{genome_assembly}/genome_seq/genomic.fna.path",genome_assembly=assembly_list)


# Ne se lance que sur un noeud (-j = 1) sinon bug
rule download_NCBI_genome:
    # Telecharge le genome du NCBI
    params:
        symlink_directory = pathGTDriftData + "genome_assembly/{genome_assembly}/genome_seq/"
    output:
        genome_path = pathGTDriftData + "genome_assembly/{genome_assembly}/genome_seq/genomic.fna.path"
    shell:
        "{pathGTDriftScripts}analyses/collecting_genome_annotation_pixi/download_genome.sh {wildcards.genome_assembly} {output.genome_path} {params.symlink_directory} {storagetype}"


# Ne se lance que sur un noeud (-j = 1) sinon bug
rule download_NCBI_protein:
    # Telecharge le fichier de sequence proteique du NCBI
    params:
        symlink_directory = pathGTDriftData + "genome_assembly/{genome_assembly}/annotation/"
    output:
        prot_path =  pathGTDriftData + "genome_assembly/{genome_assembly}/annotation/protein.faa"
    shell:
        "{pathGTDriftScripts}analyses/collecting_genome_annotation_pixi/download_protein.sh {wildcards.genome_assembly} {output.prot_path} {params.symlink_directory}"


# Ne se lance que sur un noeud (-j = 1) sinon bug
rule download_NCBI_cds:
    # Telecharge le fichier de sequence CDS du NCBI
    params:
        symlink_directory = pathGTDriftData + "genome_assembly/{genome_assembly}/annotation/"
    output:
        cds_path =  pathGTDriftData + "genome_assembly/{genome_assembly}/annotation/cds_from_genomic.fna"
    shell:
        "{pathGTDriftScripts}analyses/collecting_genome_annotation_pixi/download_cds.sh {wildcards.genome_assembly} {output.cds_path} {params.symlink_directory}"


# Ne se lance que sur un noeud (-j = 1) sinon bug
rule download_NCBI_annotation:
    # Telecharge les annotations du NCBI
    params:
        symlink_directory = pathGTDriftData + "genome_assembly/{genome_assembly}/annotation/"
    output:
        gff_path = pathGTDriftData + "genome_assembly/{genome_assembly}/annotation/genomic.gff",
        gtf_path = pathGTDriftData + "genome_assembly/{genome_assembly}/annotation/genomic.gtf"
    shell:
        "{pathGTDriftScripts}analyses/collecting_genome_annotation_pixi/download_annotation.sh {wildcards.genome_assembly} {output.gff_path} {output.gtf_path} {params.symlink_directory}"



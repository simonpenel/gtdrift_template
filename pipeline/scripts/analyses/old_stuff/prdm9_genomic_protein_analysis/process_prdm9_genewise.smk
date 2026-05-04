import os
DOMAIN = ['KRAB', 'SET', 'SSXRD', 'ZF']
EXON = ['CDS_exon2', 'CDS_exon3', 'CDS_exon4', 'CDS_exon5', 'CDS_exon6', 'CDS_exon7',
        'CDS_exon8', 'CDS_exon9', 'CDS_exon10', 'CDS_exon11']

# Function to load JSON files
def load_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

# Assign environment variables
globals().update(load_json("../environment_path.json"))

with open("data/resources/organisms_data") as reader:
    """
    Get the list of curated ones
    """
    UNCURATED = [] # Assemblies without annotation and protein sequence
    CURATED = [] # Assemblies with annotation and protein sequence
    for line in reader.readlines()[1:]:
        line_data = line.strip().split('\t')
        if line_data[-1] != 'None' and line_data[3] == 'True': # if there is an existing URL and genome is curated
                CURATED.append(line_data[2])
        elif line_data[-1] != 'None':
                UNCURATED.append(line_data[2])
if CURATED == []:
        ACCESSNB = UNCURATED
elif UNCURATED == []:
        ACCESSNB = CURATED
else:
        ACCESSNB = UNCURATED + CURATED

rule all:
    """
    Get the prdm9 stats
    """
    input:
        stats_prdm9 = expand(pathGTDriftData+ "genome_assembly/{accession}/analyses/prdm9_genomic_protein/summary_table_prdm9_{accession}.csv", accession=ACCESSNB),
        krab = pathGTDriftGlobalResults+"prdm9_genomic_protein_analysis/summarized_results/krab_data.csv",
        krabzf = pathGTDriftGlobalResults+"prdm9_genomic_protein_analysis/summarized_results/krabzf_data.csv",
        zf = pathGTDriftGlobalResults+"prdm9_genomic_protein_analysis/summarized_results/zf_count.csv",
        table = pathGTDriftGlobalResults+"prdm9_genomic_protein_analysis/summarized_results/table_prdm9.csv",
        ##hits = expand("results/{accession}/predicted_proteins/homologue_hits.faa", accession=ACCESSNB)

rule copy_results:
    """
    Copy results
    """
    input:
        #expand("results/{accession}/Result_tables/summary_table_prdm9_{accession}.csv", accession=ACCESSNB),
        "results/{accession}/Result_tables/summary_table_prdm9_{accession}.csv"
    output:
        pathGTDriftData+ "genome_assembly/{accession}/analyses/prdm9_genomic_protein/summary_table_prdm9_{accession}.csv"

    shell:
        """
        echo copy  {input} into {output};
        cp {input} {output};
        """
 

#if CURATED != []:
#        include: "module_get_faa.smk"
#        include: "module_get_gff.smk"
#include: "module_get_fna.smk" 
include: "module_prdm9_genewise.smk"


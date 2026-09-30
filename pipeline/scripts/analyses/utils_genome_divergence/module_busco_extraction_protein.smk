import json

# Function to load JSON files
def load_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

# Assign environment variables
globals().update(load_json("../environment_path.json"))

# Configuration
# -------------
configfile: "analyse.json"
configfile: "assemblies_annoted.json"

# List of assemblies
# ------------------
ACCESSNB = config["assembly_list_annoted"]

# Storage type
# -------------
# Indicates if the genome data are stored localy or on iRODs
storagetype  = config["storagetype"]


BUSCO_DIR = config["busco_dir"]
pathBUSCO = BUSCO_DIR + "protein/"

rule all:
    input:
        expand(BUSCO_DIR + "extracted_buscos/{accession}_protein_all_buscos.fa", accession = ACCESSNB)

rule clean_isoforms:
    """
    Remove isoforms from a protein fasta
    """
    priority: 1
    input:
        faa = pathGTDriftData+ "genome_assembly/{accession}/annotation/protein.faa",
        gff = pathGTDriftData+ "genome_assembly/{accession}/annotation/genomic.gff"
    output:
        clean = temp(pathGTDriftData+ "genome_assembly/{accession}/annotation/clean_protein.faa")
    shell:
        """
        python3 ../utils/python/filter_isoforms.py -f {input.faa} -g {input.gff} -o {output}
        """



rule busco_protein:
    """
    Execute BUSCO on anottated data
    """
    priority: 2
    input:
        faa = pathGTDriftData+ "genome_assembly/{accession}/annotation/clean_protein.faa"
    output:
        table = pathBUSCO + "{accession}/run_eukaryota_odb12/full_table.tsv"
    shell:
        """
        busco -i {input} -f --offline --download_path ./busco_downloads -m protein -l eukaryota_odb12 -c 1 -o {pathBUSCO}{wildcards.accession}
        """

rule extract_protein_ids:
    """
    Extract protein IDs for all BUSCO
    """
    priority: 3
    input:
        table = pathBUSCO + "{accession}/run_eukaryota_odb12/full_table.tsv"
    output:
        prots = pathBUSCO + "{accession}/extracted_protein_ids"
    shell:
        """
        python3 ../utils/python/extract_sequences_protein.py -i {input} -o {output}
        """

rule busco_extract_protein:
    """
    Extract BUSCO sequences based on protein IDs
    """
    priority: 4
    input:
        prots = pathBUSCO + "{accession}/extracted_protein_ids",
        gff = pathGTDriftData+ "genome_assembly/{accession}/annotation/genomic.gff",
        fna = pathGTDriftData+ "genome_assembly/{accession}/genome_seq/genomic.fna"
    output:
        BUSCO_DIR + "extracted_buscos/{accession}_protein_all_buscos.fa"
    shell:
        """
        mkdir -p {pathBUSCO}extracted_buscos/
        python3 ../utils/python/extract_protein_cds.py -p {input.prots} -f {input.fna} -g {input.gff} -o {pathBUSCO}extracted_buscos -a {wildcards.accession}
        cat {pathBUSCO}extracted_buscos/{wildcards.accession}_*.fasta > {output}
        """

rule get_genome_seq_fasta:
    input:
        fasta = pathGTDriftData+ "genome_assembly/{accession}/genome_seq/genomic.fna.path"
    output:
        fasta = temp(pathGTDriftData+ "genome_assembly/{accession}/genome_seq/genomic.fna")
    shell:
        """
        export  genomic=`cat {input.fasta}`
        echo "Genome sequence fasta file : $genomic"
        if [ {storagetype} == irods ];
            then
            echo "iget  /lbbeZone/home/penel/gtdrift/genome_seq/$genomic"
            ls {pathGTDriftData}"genome_assembly/{wildcards.accession}/genome_seq/"
            iget -f /lbbeZone/home/penel/gtdrift/genome_seq/$genomic {output.fasta}
        else
            ln -s {pathGTDriftData}"genome_assembly/{wildcards.accession}/genome_seq/$genomic {output.fasta}"
        fi    
        """

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
configfile: "assemblies_not_annoted.json"

# Storage type
# -------------
# Indicates if the genome data are stored localy or on iRODs
storagetype  = config["storagetype"]

# List of assemblies
# ------------------
ALL_ACCESSNB = config["assembly_list_not_annoted"]
PARTITION = config["partition"]
NB_PARTITIONS = config["nb_partitions"]

NB_ACCES = len(ALL_ACCESSNB)
SIZE_PARTITION = NB_ACCES // NB_PARTITIONS 
print({NB_ACCES})
print({NB_PARTITIONS})
print({SIZE_PARTITION})

print("debug")
if PARTITION == NB_PARTITIONS:
    ACCESSNB = ALL_ACCESSNB[(PARTITION - 1) * SIZE_PARTITION : NB_ACCES]
else :
    ACCESSNB = ALL_ACCESSNB[(PARTITION - 1) * SIZE_PARTITION : PARTITION * SIZE_PARTITION]


# Name of global results directory.
# ---------------------------------
# The directory is located in pathGTDriftGlobalResults
GLOBAL_RESULTS = config["global_analyse_dir_name"]

# Name of genome specific results directory.
# ------------------------------------------
# The directory is located in genome_assembly/{accession}/analyses/
GENOME_RESULTS = config["analyse_dir_name"]


#BUSCO_DIR = config["busco_dir"]
#pathBUSCO = BUSCO_DIR + "genomic/"




rule all:
    input:
        expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
            + "results/run_eukaryota_odb12/full_table.tsv",
            accession = ACCESSNB),
        expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
            + "single_copy_busco_sequences.gff",
            accession = ACCESSNB),
        expand(
            pathGTDriftGlobalResults + GLOBAL_RESULTS 
            + "extracted_buscos/{accession}_all_buscos.fa",
            accession = ACCESSNB)

rule busco_genomic:
    """
    Execute BUSCO on unanottated data
    """
    priority: 1
    input:
        fna = pathGTDriftData+ "genome_assembly/{accession}/genome_seq/genomic.fna"
    output:
        table = pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
            + "results/run_eukaryota_odb12/full_table.tsv"
    shell:
        """
        busco -i {input.fna} -f --offline --download_path ./busco_downloads -m genome -l eukaryota_odb12 -c 1 --out_path {pathGTDriftData}genome_assembly/{wildcards.accession}/analyses/{GENOME_RESULTS} -o results 
        #busco -i {input.fna} -f  --download_path ./busco_downloads -m genome -l eukaryota_odb12 -c 1 -o {pathGTDriftData}genome_assembly/{wildcards.accession}/analyses/{GENOME_RESULTS}
 
        """

rule concatenate_gffs_genomic:
    """
    Concatenate gffs from genomic BUSCO execution
    Checks if single or multi-copy BUSCO genes were found, and if so concatenates them.
    If none were found, touch an empty file.
    Also deletes sizeable log and temp files.
    """
    input:
        table = pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
            + "results/run_eukaryota_odb12/full_table.tsv"
    output:
        gff = pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS
            + "single_copy_busco_sequences.gff"
    shell:
        """
        if [ $(ls '{pathGTDriftData}genome_assembly/{wildcards.accession}/analyses/{GENOME_RESULTS}run_eukaryota_odb12/busco_sequences/single_copy_busco_sequences/*.gff' 2> /dev/null |wc -l) -gt 0 ]; then
            find {pathGTDriftData}genome_assembly/{wildcards.accession}/analyses/{GENOME_RESULTS}run_eukaryota_odb12/busco_sequences/single_copy_busco_sequences/*.gff -type f -print -exec cat {{}} \; > {output}
        fi
        if [ $(ls '{pathGTDriftData}genome_assembly/{wildcards.accession}/analyses/{GENOME_RESULTS}/run_eukaryota_odb12/busco_sequences/multi_copy_busco_sequences/*.gff' 2> /dev/null |wc -l) -gt 0 ]; then
            find {pathGTDriftData}genome_assembly/{wildcards.accession}/analyses/{GENOME_RESULTS}/run_eukaryota_odb12/busco_sequences/multi_copy_busco_sequences/*.gff -type f -print -exec cat {{}} \; >> {output}
            for p in $(ls {pathGTDriftData}genome_assembly/{wildcards.accession}/analyses/{GENOME_RESULTS}/run_eukaryota_odb12/busco_sequences/multi_copy_busco_sequences/*.gff);
            do
                echo $(basename $p) >> {pathGTDriftData}genome_assembly/{wildcards.accession}/analyses/{GENOME_RESULTS}/multi_copy_buscos
            done
        fi
        if [ $(ls '{pathGTDriftData}genome_assembly/{wildcards.accession}/analyses/{GENOME_RESULTS}/run_eukaryota_odb12/busco_sequences/single_copy_busco_sequences/*.gff' 2> /dev/null|wc -l) -gt 0 ] || [ $(ls '{pathGTDriftData}genome_assembly/{wildcards.accession}/analyses/{GENOME_RESULTS}/run_eukaryota_odb12/busco_sequences/multi_copy_busco_sequences/*.gff' 2> /dev/null |wc -l) -gt 0 ]; then
            echo "GFFs concatenated."
        else
            echo "No BUSCO found for {wildcards.accession}"
            touch {output}
        fi
 #       if [ -f '{pathGTDriftData}genome_assembly/{wildcards.accession}/analyses/{GENOME_RESULTS}/run_eukaryota_odb12/miniprot_output/ref.mpi' ]; then
 #           rm {pathGTDriftData}genome_assembly/{wildcards.accession}/analyses/{GENOME_RESULTS}/run_eukaryota_odb12/miniprot_output/ref.mpi
 #       fi
 #       if [ -f '{pathGTDriftData}genome_assembly/{wildcards.accession}/analyses/{GENOME_RESULTS}/logs/miniprot_align_eukaryota_odb12_out.log' ]; then
 #           rm {pathGTDriftData}genome_assembly/{wildcards.accession}/analyses/{GENOME_RESULTS}/logs/miniprot_align_eukaryota_odb12_out.log
 #       fi
 #       if [ -d '{pathGTDriftData}genome_assembly/{wildcards.accession}/analyses/{GENOME_RESULTS}/tmp/' ]; then
 #           rm -r {pathGTDriftData}genome_assembly/{wildcards.accession}/analyses/{GENOME_RESULTS}/tmp/
 #       fi
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

rule busco_extract_genomic:
    """
    Extract BUSCO sequences from the genomic fasta
    """
    priority: 3
    input:
        fna = pathGTDriftData + "genome_assembly/{accession}/genome_seq/genomic.fna",
        gff = pathGTDriftData 
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
        + "single_copy_busco_sequences.gff"
    output:
        pathGTDriftGlobalResults + GLOBAL_RESULTS 
        + "extracted_buscos/{accession}_all_buscos.fa"
    shell:
        """
        python3 python/extract_genomic_cds.py -f {input.fna} -g {input.gff} -o {pathGTDriftData}genome_assembly/{wildcards.accession}/analyses/{GENOME_RESULTS}extracted_buscos -a {wildcards.accession}
        if compgen -G '{pathGTDriftData}genome_assembly/{wildcards.accession}/analyses/{GENOME_RESULTS}extracted_buscos/{wildcards.accession}_*.fasta' > /dev/null; then
            cat {pathGTDriftData}genome_assembly/{wildcards.accession}/analyses/{GENOME_RESULTS}extracted_buscos/{wildcards.accession}_*.fasta > {output}
        else
            touch {output}
        fi
        """        
## Module for prdm9 analysis on protein data
## Date : Decembre 2024
## Authors :


# Configuration
# -------------


configfile: "config.json"


ACCESSNB = config["assembly_list"]
DOMAIN = ["KRAB", "SET", "SSXRD", "ZF"]


# Function to load JSON files
# ---------------------------
def load_json(file_path):
    with open(file_path, "r") as file:
        return json.load(file)


# Assign environment variables
# ----------------------------
globals().update(load_json("../environment_path.json"))



# Rules
# -----

# -------------------------------------------------------------------
# summary
# concatenate all the candidate search in human PRDM family summaries
# into a single file.
# -------------------------------------------------------------------
#rule summary:
#    """
#    Concatenation of each proteome blastp results.
#    """
#    input:
#        # summaries the search in human PRDM family for each assemblies
#        expand(
#            pathGTDriftData
#            + "genome_assembly/{accession}/analyses/prdm9_prot/blastp.txt",
#            accession=ACCESSNB,
#        ),
#    output:
#        # concatenation of all summaries
#        pathGTDriftGlobalResults
#        + "analyses_summaries/BLASTP_results/blastp_summary.txt",
#    shell:
#        """
#        cat {input} > {output}
#        """


# --------------------------------------------------------------
# blastp_results
# create a table of the different candidate for all  assemblies.
# ;Taxid;Species_name;Protein ID;Bit score;Ratio;SET;KRAB;SSXRD;ZF;Best_Match
# Warning : the python script use  the file ../data_results_per_assembly/organisms_data
# --------------------------------------------------------------
#rule blastp_results:
#    """
#    Writing a table from the concatenation
#    """
#    input:
#        # concatenation of all summaries
#        pathGTDriftGlobalResults
#        + "analyses_summaries/BLASTP_results/blastp_summary.txt",
#    output:
#        # table of the different candidate for all  assemblies
#        pathGTDriftGlobalResults
#        + "analyses_summaries/BLASTP_results/blastp_results.csv",
#    shell:
#        (
#            "python3 "
#            + pathGTDriftScripts
#            + "analyses/prdm9_protein_analysis//python/blastp_table.py -i "
#            + pathGTDriftGlobalResults
#        )


# -----------------------------------------------------------------
# taxonomy
# Creation of a table associating a genome accession number to its
# complete taxonomy.
# Warning : the python script use  the file ../data_results_per_assembly/organisms_data
# -----------------------------------------------------------------
#rule taxonomy:
#    """
#    Creation of a table associating a genome accession number to its complete taxonomy
#    """
#    input:
#        # concatenation of all summaries
#        pathGTDriftGlobalResults
#        + "analyses_summaries/BLASTP_results/blastp_summary.txt",
#    output:
#        # full taxonomy for each assemblies
#        pathGTDriftGlobalResults + "sorted_taxonomy.csv",
#    shell:
#        (
#            "python3 "
#            + pathGTDriftScripts
#            + "analyses/prdm9_protein_analysis/python/taxonomy.py -i "
#            + pathGTDriftData
#            + " -o "
#            + pathGTDriftGlobalResults
#        )

# -----------------------------------------------------------------
# create_global_krab_table
# write a summary of krab domain in the proteome.
# Output format:
# Accession;Protein List;KRAB nb
# -----------------------------------------------------------------
rule create_global_krab_table:
    """
    Creation of global krab result table 
    """
    input:
        krab_tabulated = expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/tbl/KRAB_tabulated",
            accession=ACCESSNB,
        ),
    output:
        pathGTDriftGlobalResults + "analyses_summaries/table_results/krab_data.csv",
    shell:
        (
            "python3 "
            + pathGTDriftScripts
            + "/analyses/prdm9_protein_analysis/python/krab.py  -i "
            + pathGTDriftData
            + " -o  {output}"
        )
        
# -----------------------------------------------------------------
# create_global_krabzf_table
# write a summary of krab domain in the proteome.
# Output format:
# Accession;KRAB+ZF protein list;KRAB+ZF nb
# -----------------------------------------------------------------
rule create_global_krabzf_table:
    """
    Creation of global krab and zf result table
    """
    input:
        krab = pathGTDriftGlobalResults + "analyses_summaries/table_results/krab_data.csv",
        zf_tabulated = expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/tbl/ZF_tabulated",
            accession=ACCESSNB,
        ),
    output:
        pathGTDriftGlobalResults + "analyses_summaries/table_results/krabzf_data.csv",
    shell:
        (
            "python3 "
            + pathGTDriftScripts
            + "/analyses/prdm9_protein_analysis/python/krabzf.py  -i "
            + pathGTDriftData
            + " -k  {input.krab}"
            + " -o  {output}"
        )      
        
# -----------------------------------------------------------------
# create_global_zf_table
# write a summary of zf protein with more than 5 zf domains.
# Output format:
# Accession;5+ ZF
# -----------------------------------------------------------------
rule create_global_zf_table:
    """
    Creation of global zf count table
    """
    input:
        zf_domains = expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/ZF_domains_summary",
            accession=ACCESSNB,
        ),
    output:
        pathGTDriftGlobalResults + "analyses_summaries/table_results/zf_count.csv",
    shell:
        (
            "python3 "
            + pathGTDriftScripts
            + "/analyses/prdm9_protein_analysis/python/zf_analysis.py  -i "
            + pathGTDriftData
            + " -o  {output}"
        )             


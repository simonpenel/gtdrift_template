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


if config["mode"] == "guix":
    RUNCMD = "guix shell hmmer -- "
else:
    RUNCMD = ""


# Rules
# -----


# get_blastdb: Build a blast database from the proteins wich should 
# be stored in the genome_assembly directory. These data were 
# collected with the "collecting_genome_annotation" pipeline.
# Use of formatdb instead makeblastdb because last makeblastdb 
# version doestn work on beegfs. Alternatively, it is possible to use
# an old version of makeblasdb by modify the PATH in your bash 
# environment.
# -------------------------------------------------------------------
rule get_blastdb:
    """
    Genere a blast db
    WARNING: formatdb is used and not makeblastdb. Suffixes and options are different with makeblastdb
    """
    input:
        # the protein fasta file
        fasta=pathGTDriftData + "genome_assembly/{accession}/annotation/protein.faa",
    output:
        # some databse index 
        # (we dont need to specify not all of them, TODO : use a directory instead)
        psq=pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.psq",
        pin=pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.pin",
        phr=pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.phr",
    shell:
        #"makeblastdb -in {input.fasta} -title protdb -out "+pathGTDriftData+"genome_assembly/{wildcards.accession}/analyses/prdm9_prot/protdb -dbtype prot "
        (
            "formatdb -i {input.fasta} -t protdb -n "
            + pathGTDriftData
            + "genome_assembly/{wildcards.accession}/analyses/prdm9_prot/protdb -p T -o T"
        )


# hmm_build
# build the HMM from the Prdm9 metazoa reference alignments for each
# domain.
# ------------------------------------------------------------------
rule hmm_build:
    """
    HMM creation.
    """
    input:
        # the domain alignment in fasta format
        pathGTDriftResource
        + "ref_align/Prdm9_Metazoa_Reference_alignment/Domain_{domain}_ReferenceAlignment.fa",
    output:
        # the hmm profile for the domain
        pathGTDriftResource + "hmm_build/{domain}.hmm",
    shell:
        "{RUNCMD} hmmbuild {output} {input}"


# hmm_search
# search the  HMM profile of the domain in the proteome.
# ------------------------------------------------------
rule hmm_search:
    """
    Proteome search using the HMMs.
    """
    input:
        # the hmm profile of the domain
        model=pathGTDriftResource + "hmm_build/{domain}.hmm",
        # the proteome
        protein=pathGTDriftData + "genome_assembly/{accession}/annotation/protein.faa",
    output:
        # table of of per-sequence hits
        table=pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/tbl/{domain}",
        # table of per-domain hits
        domains=pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains",
    shell:
        "{RUNCMD} hmmsearch -E 1E-3 --domE 1E-3 --tblout {output.table} --domtblout {output.domains} --noali {input.model} {input.protein}"


# tbl_processing
# write per-sequence hits in tabular format.
# ------------------------------------------
rule tbl_processing:
    """
    Result file processing for a later use.
    """
    input:
        # per-sequence hits from hmmsearch
        pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/tbl/{domain}",
    output:
        # per-sequence hits in tabular format
        pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/tbl/{domain}_processed",
    shell:
        (
            "python3 "
            + pathGTDriftScripts
            + "analyses/prdm9_protein_analysis/python/hmmsearch_parser.py -i {input} -o {output}"
        )


# domain_processing
# write per-domain hits in tabular format.
# ----------------------------------------
rule domain_processing:
    """
    Result file processing for a later use.
    """
    input:
        # per-sequence hits in tabular format
        pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/tbl/{domain}_processed",
        # per-domain hits from hmmsearch
        domain_data=pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains",
    output:
        # per-domain hits in tabular format
        processed=pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains_processed",
        # per-domain hits in tabular format in which overlapping zinc finger domains are
        # merged to create one big domain with multiple repetitions 
        summary=pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains_summary",
    shell:
        (
            "python3  "
            + pathGTDriftScripts
            + "analyses/prdm9_protein_analysis/python/domain_parser.py -i {input.domain_data} -o {output.processed} -s {output.summary}"
        )


# function returning the path of all per-domain hits in tabular 
# format with overlapping zinc finger domains.
# -------------------------------------------------------------
def domain_done(wildcards):
    return expand(
        pathGTDriftData
        + "genome_assembly/"
        + wildcards.accession
        + "/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains_summary",
        domain=DOMAIN,
    )


# table_editing
# write prdm9 protein statistics for each assembly.
# warning : summary_table_prdm9_{accession}.csv should not be 
# confused with summary_table_{accession}.csv
# Output format:
# ;SeqID;SET Query;SET E-value;SET Score;Nb SET domains;SET domain start;SET domain end;KRAB Query;KRAB E-value;KRAB Score;Nb KRAB domains;KRAB domain start;KRAB domain end;SSXRD Query;SSXRD E-value;SSXRD Score;Nb SSXRD domains;SSXRD domain start;SSXRD domain end;ZF Query;ZF E-value;ZF Score;Nb ZF domains;ZF domain start;ZF domain end;Taxid
# ------------------------------------------------------------
rule table_editing:
    """
    Creation of a summary table of hmm_search results.
    """
    input:
        # path of all per-domain hits in tabular format with overlapping zinc finger domains
        domain_done,
    output:
        # prdm9 protein statistics for each assembly.
        # warning : summary_table_prdm9_{accession}.csv should
        # not be confused with summary_table_{accession}.csv
        pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/summary_table_prdm9_{accession}.csv",
    shell:
        (
            "python3 "
            + pathGTDriftScripts
            + "analyses/prdm9_protein_analysis/python/table_builder.py -i "
            + pathGTDriftData
            + "genome_assembly -a {wildcards.accession} -o {output}"
        )


# read_table
# Extracts the sequence selected by hmm search for an organism and
# runs a blastp analysis against the Human PRDM genes family. If the
# best match is PRDM9, the value is saved and compared to the next 
# best non-PRDM9 match. The ouput blastp.txt file contains the taxid,
# the best PRDM match,the presence/absence data for every proteic
# domain, the bit score of the blastp if the best match is PRDM9 and
# the ratio with the second best non-PRDM9 match.
# The summary_table_{accession}.csv is more detailed :
# ;Unnamed: 0;SeqID;SET Query;SET E-value;SET Score;Nb SET domains;SET domain start;SET domain end;KRAB Query;KRAB E-value;KRAB Score;Nb KRAB domains;KRAB domain start;KRAB domain end;SSXRD Query;SSXRD E-value;SSXRD Score;Nb SSXRD domains;SSXRD domain start;SSXRD domain end;ZF Query;ZF E-value;ZF Score;Nb ZF domains;ZF domain start;ZF domain end;Taxid;Best Match;Bit Score;Score ratio
# -------------------------------------------------------------------
rule read_table:
    """
    Reads each summary table and runs a blastp analysis on every candidate
    """
    input:
        # prdm9 protein statistics for each assembly.
        pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/summary_table_prdm9_{accession}.csv",
        # some of the blast database index (TODO use a directory instead)        
        psq=pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.psq",
        # some of the blast database index (TODO use a directory instead)
        pin=pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.pin",
        # some of the blast database index (TODO use a directory instead)        
        phr=pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.phr",
    output:
        # a summary of the blastp results of the search in the human PRDM family    
        pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/blastp.txt",
        # detailed results inluding the blastp results        
        pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/summary_table_{accession}.csv",
    shell:
        (
            "python3 "
            + pathGTDriftScripts
            + "analyses/prdm9_protein_analysis/python/blastp_analysis.py "
            + pathGTDriftData
            + "genome_assembly/{wildcards.accession}/analyses/prdm9_prot/summary_table_prdm9_{wildcards.accession}.csv {wildcards.accession} "
            + pathGTDriftData
            + "genome_assembly/"
        )


# summary
# concatenate all the candidate search in human PRDM family summaries
# into a unique file.
# -------------------------------------------------------------------
rule summary:
    """
    Concatenation of each proteome blastp results.
    """
    input:
        # summaries the search in human PRDM family for each assemblies
        expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/prdm9_prot/blastp.txt",
            accession=ACCESSNB,
        ),
    output:
        # concatenation of all summaries
        pathGTDriftGlobalResults
        + "analyses_summaries/BLASTP_results/blastp_summary.txt",
    shell:
        """
        cat {input} > {output}
        """


# blastp_results
# create a table of the different candidate for all  assemblies.
# ;Taxid;Species_name;Protein ID;Bit score;Ratio;SET;KRAB;SSXRD;ZF;Best_Match
# Warning : the python script use  the file ../data_results_per_assembly/organisms_data
# --------------------------------------------------------------
rule blastp_results:
    """
    Writing a table from the concatenation
    """
    input:
        # concatenation of all summaries
        pathGTDriftGlobalResults
        + "analyses_summaries/BLASTP_results/blastp_summary.txt",
    output:
        # table of the different candidate for all  assemblies
        pathGTDriftGlobalResults
        + "analyses_summaries/BLASTP_results/blastp_results.csv",
    shell:
        (
            "python3 "
            + pathGTDriftScripts
            + "analyses/prdm9_protein_analysis//python/blastp_table.py -i "
            + pathGTDriftGlobalResults
        )


# taxonomy
# Creation of a table associating a genome accession number to its
# complete taxonomy.
# Warning : the python script use  the file ../data_results_per_assembly/organisms_data
# -----------------------------------------------------------------
rule taxonomy:
    """
    Creation of a table associating a genome accession number to its complete taxonomy
    """
    input:
        # concatenation of all summaries
        pathGTDriftGlobalResults
        + "analyses_summaries/BLASTP_results/blastp_summary.txt",
    output:
        # full taxonomy for each assemblies
        pathGTDriftGlobalResults + "sorted_taxonomy.csv",
    shell:
        (
            "python3 "
            + pathGTDriftScripts
            + "analyses/prdm9_protein_analysis/python/taxonomy.py -i "
            + pathGTDriftData
            + " -o "
            + pathGTDriftGlobalResults
        )


# create_table
# Creation of multiple result table using blastp results and hmm 
# search results.
# TODO: sorted_taxonomy does not seem to be used in the scripts.
# What is the purpose of this input?
# ---------------------------------------------------------------
rule create_table:
    """
    Creation of multiple result table using blastp results and hmm search results
    """
    input:
        pathGTDriftGlobalResults
        + "analyses_summaries/BLASTP_results/blastp_results.csv",
        pathGTDriftGlobalResults + "sorted_taxonomy.csv",
    output:
        pathGTDriftGlobalResults + "analyses_summaries/table_results/krab_data.csv",
        pathGTDriftGlobalResults + "analyses_summaries/table_results/krabzf_data.csv",
        pathGTDriftGlobalResults + "analyses_summaries/table_results/zf_count.csv",
    shell:
        (
            "python3 "
            + pathGTDriftScripts
            + "/analyses/prdm9_protein_analysis/python/krab.py -i "
            + pathGTDriftData
            + " -o "
            + pathGTDriftGlobalResults
            + "\
                                                                                                                                                && python3 "
            + pathGTDriftScripts
            + "/analyses/prdm9_protein_analysis/python/krabzf.py  -i "
            + pathGTDriftData
            + " -o "
            + pathGTDriftGlobalResults
            + "\
                                                                                                                                                && python3 "
            + pathGTDriftScripts
            + "analyses/prdm9_protein_analysis/python/zf_analysis.py  -i "
            + pathGTDriftData
            + " -o "
            + pathGTDriftGlobalResults
        )

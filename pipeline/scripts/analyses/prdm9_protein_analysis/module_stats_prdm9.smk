## Module for prdm9 analysis on protein data
## Date : Decembre 2024
## Authors : L. Duret 
## Purpose: Search for PRDM9 sequences in the proteome.

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


# -------------------------------------------------------------------
# get_blastdb
# Build a blast database from the proteins wich should be stored in
# the genome_assembly directory. These data were collected with the
# "collecting_genome_annotation" pipeline. Use of formatdb instead
# makeblastdb because last makeblastdb version doestn work on beegfs.
# Alternatively, it is possible to use an old version of makeblasdb
# by modify the PATH in your bash  environment.
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


# ------------------------------------------------------------------
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


# ------------------------------------------------------
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


# ------------------------------------------
# formating_hmm_sequence_hit
# write per-sequence hits in tabular format.
# ------------------------------------------
rule formating_hmm_sequence_hit:
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
        + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/tbl/{domain}_tabulated",
    script:
        "python/hmmsearch_parser.py"

# ----------------------------------------
# formating_hmm_domain_hit
# write per-domain hits in tabular format.
# ----------------------------------------
rule formating_hmm_domain_hit:
    """
    Result file processing for a later use.
    """
    input:
        # per-domain hits from hmmsearch
        per_domain=pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains",
    output:
        # per-domain hits in tabular format
        tabulated_per_domain=pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains_tabulated",
        # per-domain hits in tabular format in which overlapping zinc finger domains are
        # merged to create one big domain with multiple repetitions 
        domain_summary=pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains_summary",
    script:
        "python/domain_parser.py"
#    shell:
#        (
#            "python3  "
#            + pathGTDriftScripts
#            + "analyses/prdm9_protein_analysis/python/domain_parser.py -i {input.domain_data} -o #{output.processed} -s {output.summary}"
#        )


# -------------------------------------------------------------
# function returning the path of all per-domain hits in tabular
# format with overlapping zinc finger domains.
# -------------------------------------------------------------
#def domain_done(wildcards):
#    return expand(
#        pathGTDriftData
#        + "genome_assembly/"
#        + wildcards.accession
#        + "/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains_summary",
#        domain=DOMAIN,
#    )

def accession_nb(wildcards):
    return wildcards.accession

# ------------------------------------------------------------
# summarize_hmm_results
# write results of hmm search  for each assembly.
# warning : summary_hmmsearch_prdm9_{accession}.csv should not be
# confused with summary_hmmsearch_prdm9_with_paralog_check_{accession}.csv
# Output format:
# ;SeqID;SET Query;SET E-value;SET Score;Nb SET domains;SET domain start;SET domain end;KRAB Query;KRAB E-value;KRAB Score;Nb KRAB domains;KRAB domain start;KRAB domain end;SSXRD Query;SSXRD E-value;SSXRD Score;Nb SSXRD domains;SSXRD domain start;SSXRD domain end;ZF Query;ZF E-value;ZF Score;Nb ZF domains;ZF domain start;ZF domain end;Taxid
# ------------------------------------------------------------
rule summarize_hmm_results:
    """
    Creation of a summary table of hmm_search results.
    """
    input:
        # path of all per-domain hits in tabular format with overlapping zinc finger domains
        #accession=accession_nb,
        # tabulated results of hmm search on SET domain
        SET_per_sequence_tabulated=pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/tbl/SET_tabulated",
        SET_per_domain_tabulated=pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/SET_domains_tabulated",
        KRAB_per_domain_tabulated=pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/KRAB_domains_tabulated",
        ZF_per_domain_tabulated=pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/ZF_domains_tabulated", 
        SSXRD_per_domain_tabulated=pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/SSXRD_domains_tabulated",               
        SET_per_domain_summary=pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/SET_domains_summary",
        KRAB_per_domain_summary=pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/KRAB_domains_summary",
        ZF_per_domain_summary=pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/ZF_domains_summary", 
        SSXRD_per_domain_summary=pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/SSXRD_domains_summary",          

    output:
        # prdm9 protein statistics for each assembly.
        # warning : summary_hmmsearch_prdm9_{accession}.csv should
        # not be confused with summary_hmmsearch_prdm9_with_paralog_check_{accession}.csv
        pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/summary_hmmsearch_prdm9_{accession}.csv",
    params:    
        accession=accession_nb,  
    script:
        "python/table_builder.py"
#    shell:
#        (
#            "python3 "
#            + pathGTDriftScripts
#            + "analyses/prdm9_protein_analysis/python/table_builder.py -i "
#            + pathGTDriftData
#            + "genome_assembly -a {wildcards.accession} -o {output}"
#        )


# -------------------------------------------------------------------
# prdm_paralog_check
# Extracts the sequence selected by hmm search for an organism and
# runs a blastp analysis against the Human PRDM genes family. If the
# best match is PRDM9, the value is saved and compared to the next
# best non-PRDM9 match. The ouput blastp.txt file contains the taxid,
# the best PRDM match,the presence/absence data for every proteic
# domain, the bit score of the blastp if the best match is PRDM9 and
# the ratio with the second best non-PRDM9 match.
# The summary_hmmsearch_prdm9_with_paralog_check_{accession}.csv is more detailed :
# ;Unnamed: 0;SeqID;SET Query;SET E-value;SET Score;Nb SET domains;SET domain start;SET domain end;KRAB Query;KRAB E-value;KRAB Score;Nb KRAB domains;KRAB domain start;KRAB domain end;SSXRD Query;SSXRD E-value;SSXRD Score;Nb SSXRD domains;SSXRD domain start;SSXRD domain end;ZF Query;ZF E-value;ZF Score;Nb ZF domains;ZF domain start;ZF domain end;Taxid;Best Match;Bit Score;Score ratio
# -------------------------------------------------------------------
rule prdm_paralog_check:
    """
    Reads each summary table and runs a blastp analysis on every candidate
    """
    input:
        # prdm9 protein statistics for each assembly.
        pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/summary_hmmsearch_prdm9_{accession}.csv",
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
        blast=pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/blastp.txt",
        # detailed results inluding the blastp results        
        table=pathGTDriftData
        + "genome_assembly/{accession}/analyses/prdm9_prot/summary_hmmsearch_prdm9_with_paralog_check_{accession}.csv",
    shell:
        (
            "python3 "
            + pathGTDriftScripts
            + "analyses/prdm9_protein_analysis/python/blastp_analysis.py "
            + pathGTDriftData
            + "genome_assembly/{wildcards.accession}/analyses/prdm9_prot/summary_hmmsearch_prdm9_{wildcards.accession}.csv {wildcards.accession} "
            + pathGTDriftData
            + "genome_assembly/  {output.table}"
        )



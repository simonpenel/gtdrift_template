## Module for prdm9 analysis on protein data
## Date : Decembre 2024
## Authors : L. Duret, A. Raimbault, S. Penel 
## Purpose: Search for PRDM9 sequences in the proteome.

# Configuration
# -------------


configfile: "config.json"

if config["mode"] == "guix":
    RUNCMD = "guix shell hmmer -- "
else:
    RUNCMD = ""

DOMAIN_HMM_DIR  = config["domain_hmm_dir"]

#ACCESSNB = config["assembly_list"]
#DOMAIN = ["KRAB", "SET", "SSXRD", "ZF"]


# Function to load JSON files
# ---------------------------
def load_json(file_path):
    with open(file_path, "r") as file:
        return json.load(file)


# Assign environment variables
# ----------------------------
#globals().update(load_json("../environment_path.json"))




# Rules
# -----


# -------------------------------------------------------------------
# build_proteome_blastdb
# Build a blast database from the proteome wich should be stored in
# the genome_assembly directory. These data were collected with the
# "collecting_genome_annotation" pipeline. Use of formatdb instead
# makeblastdb because last makeblastdb version doestn work on beegfs.
# Alternatively, it is possible to use an old version of makeblasdb
# by modify the PATH in your bash  environment.
# -------------------------------------------------------------------
#rule build_proteome_blastdb:
#    """
#    Genere a blast db 
#    """
#    input:
#        # the protein fasta file
#        fasta=pathGTDriftData + "genome_assembly/{accession}/annotation/protein.faa",
#    output:
#        # the proteome database
#        dbprot=directory(pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/dbprot"),
#    shell:
#        "mkdir {output.dbprot} && formatdb -i {input.fasta} -t protdb -n {output.dbprot}/protdb -p T -o T"


# ------------------------------------------------------------------
# hmm_build
# build the HMM from the reference alignments for each
# domain.
# ------------------------------------------------------------------
rule hmm_build:
    """
    HMM creation.
    """
    input:
        # the domain alignment in fasta format
        pathGTDriftResource
        + "ref_align/" + DOMAIN_HMM_DIR + "{domain}.fa",
    output:
        # the hmm profile for the domain
        pathGTDriftResource + "hmm_build/" + DOMAIN_HMM_DIR + "{domain}.hmm",
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
        model=pathGTDriftResource + "hmm_build/" + DOMAIN_HMM_DIR + "{domain}.hmm",
        # the proteome
        protein=pathGTDriftData + "genome_assembly/{accession}/annotation/protein.faa",
    output:
        # table of of per-sequence hits
        table=pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/tbl/{domain}",
        # table of per-domain hits
        domains=pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "/hmm_search/domtbl/{domain}_domains",
    shell:
        "{RUNCMD} hmmsearch -E 1E-3 --domE 1E-3 --tblout {output.table} --domtblout {output.domains} --noali {input.model} {input.protein}"


# ------------------------------------------
# formating_hmm_sequence_hit
# write per-sequence hits in tabular format.
# ------------------------------------------
rule formating_hmm_sequence_hit:
    """
    Formating per-sequence hits from hmmsearch.
    """
    input:
        # per-sequence hits from hmmsearch
        pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/tbl/{domain}",
    output:
        # per-sequence hits in tabular format
        pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/tbl/{domain}_tabulated",
    script:
        "python/hmmsearch_parser.py"

# ----------------------------------------
# formating_hmm_domain_hit
# write per-domain hits in tabular format.
# ----------------------------------------
rule formating_hmm_domain_hit:
    """
    Formating per-domain hits from hmmsearch.
    """
    input:
        # per-domain hits from hmmsearch
        per_domain=pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/domtbl/{domain}_domains",
    output:
        # per-domain hits in tabular format
        #tabulated_per_domain=temp(pathGTDriftData
        #+ "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains_tabulated"),
        # per-domain hits in tabular format in which overlapping zinc finger domains are
        # merged to create one big domain with multiple repetitions 
        # (ZF_domain_summary will differ from ZF_domain_tabulated,
        # but there will no difference between the 2 ouputs for the
        # other domains)
        domain_summary=pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/domtbl/{domain}_domains_summary",
    script:
        "python/domain_parser.py"

# ------------------------------------------------------------
# Function sending the accession number
# ------------------------------------------------------------
def accession_nb(wildcards):
    return wildcards.accession

# ------------------------------------------------------------
# summarize_hmm_results
# write results of hmm search  for each assembly.
# warning : summary_hmmsearch_prdm9_{accession}.csv should not be
# confused with summary_hmmsearch_prdm9_with_paralog_check_{accession}.csv
# Output format:
# ;SeqID;SET Query;SET E-value;SET Score;Nb SET domains;SET domain start;SET domain end;KRAB Query;KRAB E-value;KRAB Score;Nb KRAB domains;KRAB domain start;KRAB domain end;SSXRD Query;SSXRD E-value;SSXRD Score;Nb SSXRD domains;SSXRD domain start;SSXRD domain end;ZF Query;ZF E-value;ZF Score;Nb ZF domains;ZF domain start;ZF domain end;Taxid
#### A AMELIORER : LORSQU'UN DOMAINE A PLUSIEUR HITS, LE DERNIER EST RETENU ALORS QUE C'EST LE MOINS STATISTIQUEMENT SIGNIFIANT
# ------------------------------------------------------------
#rule summarize_hmm_results:
#    """
#    Creation of a summary table of hmm_search results.
#    """
#    input:
#        # organisms_data file
#        organisms_file=pathGTDriftData + "organisms_data",
#        # tabulated results of hmm search on SET domain
#        SET_per_sequence_tabulated=pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/tbl/SET_tabulated",
#        # path of all per-domain hits in tabular format with overlapping zinc finger domains                     
#        SET_per_domain_summary=pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/SET_domains_summary",
#        KRAB_per_domain_summary=pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/KRAB_domains_summary",
#        ZF_per_domain_summary=pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/ZF_domains_summary", 
#        SSXRD_per_domain_summary=pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/SSXRD_domains_summary",          

#    output:
#        # prdm9 protein statistics for each assembly.
#        # warning : summary_hmmsearch_prdm9_{accession}.csv should
#        # not be confused with summary_hmmsearch_prdm9_with_paralog_check_{accession}.csv
#        pathGTDriftData
#        + "genome_assembly/{accession}/analyses/prdm9_prot/summary_hmmsearch_prdm9_{accession}.csv",
#    params:    
#        accession=accession_nb,  
#    script:
#        "python/table_builder.py"



configfile: "config.json"
ACCESSNB  = config["assembly_list"]
DOMAIN = ['KRAB', 'SET', 'SSXRD', 'ZF']


# Function to load JSON files
def load_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

# Assign environment variables
globals().update(load_json("../environment_path.json"))


if config["mode"] == "guix":
    RUNCMD="guix shell hmmer -- "
else:
    RUNCMD=""

rule get_blastdb:
    """
    Genere a blast db
    WARNING: formatdb is used and not makeblastdb. Suffixes and options are different with makeblastdb
    """
    input:
        fasta = pathGTDriftData+ "genome_assembly/{accession}/annotation/protein.faa"
    output:    
        psq= pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.psq",
    #    psi= pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.psi",
    #    psd= pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.psd",
        pin= pathGTDriftData+ "genome_assembly/{accession}/analyses/prdm9_prot/protdb.pin",
        phr= pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.phr"        
    shell:
       #"makeblastdb -in {input.fasta} -title protdb -out "+pathGTDriftData+"genome_assembly/{wildcards.accession}/analyses/prdm9_prot/protdb -dbtype prot "
       "formatdb -i {input.fasta} -t protdb -n "+pathGTDriftData+"genome_assembly/{wildcards.accession}/analyses/prdm9_prot/protdb -p T -o T"

         
rule hmm_build:
    """
    HMM creation.
    """
    input:
        pathGTDriftResource + "ref_align/Prdm9_Metazoa_Reference_alignment/Domain_{domain}_ReferenceAlignment.fa"
    output:
       	pathGTDriftResource + "hmm_build/{domain}.hmm"
    shell:
        "{RUNCMD} hmmbuild {output} {input}"
 
rule hmm_search:
    """
    Proteome search using the HMMs.
    """
    input:
        model=pathGTDriftResource + "hmm_build/{domain}.hmm",
        protein=pathGTDriftData + "genome_assembly/{accession}/annotation/protein.faa"
    output:
        table = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/tbl/{domain}",
        domains = pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains"
    shell:
        "{RUNCMD} hmmsearch -E 1E-3 --domE 1E-3 --tblout {output.table} --domtblout {output.domains} --noali {input.model} {input.protein}"


rule tbl_processing:
    """
    Result file processing for a later use.
    """
    input:
        pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/tbl/{domain}"
    output:
        pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/tbl/{domain}_processed"
    shell:
        "python3 "+pathGTDriftScripts+"analyses/prdm9_protein_analysis/python/hmmsearch_parser.py -i {input} -o {output}"

rule domain_processing:
    """
    Result file processing for a later use.
    """
    input:
        pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/tbl/{domain}_processed",
        domain_data=pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains"
    output:
        processed=pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains_processed",
        summary=pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains_summary"
    shell:
       "python3  "+pathGTDriftScripts+"analyses/prdm9_protein_analysis/python/domain_parser.py -i {input.domain_data} -o {output.processed} -s {output.summary}"


def domain_done(wildcards):
    return expand(pathGTDriftData + "genome_assembly/"+ wildcards.accession + "/analyses/prdm9_prot/hmm_search/domtbl/{domain}_domains_summary" , domain=DOMAIN)
    #return expand("results/" + wildcards.accession + "/hmm_search/domtbl/{domain}_domains_summary" , domain=DOMAIN)

rule table_editing:
    """
    Creation of a summary table of hmm_search results.
    """
    input:
        domain_done,
    output:
        pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/summary_table_prdm9_{accession}.csv"
    shell:
        "python3 "+pathGTDriftScripts+"analyses/prdm9_protein_analysis/python/table_builder.py -i "+ pathGTDriftData + "genome_assembly -a {wildcards.accession} -o {output}"


rule read_table:
    """
    Reads each summary table and runs a blastp analysis on every candidate
    """
    input:
        pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/summary_table_prdm9_{accession}.csv",
        psq= pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.psq",
        #psi= pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.psi",
        #psd= pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.psd",
        pin= pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.pin",
        phr= pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/protdb.phr"   
        
    output:
        pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/blastp.txt",
        pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/summary_table_{accession}.csv"
    shell:
        "python3 "+pathGTDriftScripts+"analyses/prdm9_protein_analysis/python/blastp_analysis.py "+ pathGTDriftData + "genome_assembly/{wildcards.accession}/analyses/prdm9_prot/summary_table_prdm9_{wildcards.accession}.csv {wildcards.accession} "+ pathGTDriftData + "genome_assembly/"

rule summary:
    """
    Concatenation of each proteome blastp results.
    """
    input: 
        expand(pathGTDriftData + "genome_assembly/{accession}/analyses/prdm9_prot/blastp.txt", accession=ACCESSNB)
    output:
        pathGTDriftGlobalResults + "analyses_summaries/BLASTP_results/blastp_summary.txt"
    shell:
        """
        cat {input} > {output}
        """

rule blastp_results:
    """
    Writing a table from the concatenation
    """
    input:
        pathGTDriftGlobalResults + "analyses_summaries/BLASTP_results/blastp_summary.txt"
    output:
        pathGTDriftGlobalResults + "analyses_summaries/BLASTP_results/blastp_results.csv",
    shell:
        "python3 "+pathGTDriftScripts+"analyses/prdm9_protein_analysis//python/blastp_table.py -i " + pathGTDriftGlobalResults

rule taxonomy:
    """
    Creation of a table associating a genome accession number to its complete taxonomy
    """
    input:
        pathGTDriftGlobalResults + "analyses_summaries/BLASTP_results/blastp_summary.txt"
    output:
        pathGTDriftGlobalResults + "sorted_taxonomy.csv"
    shell:
        "python3 "+pathGTDriftScripts+"analyses/prdm9_protein_analysis/python/taxonomy.py -i " +  pathGTDriftData + " -o " +  pathGTDriftGlobalResults

rule create_table:
    """
    Creation of multiple result table using blastp results and hmm search results
    """
    input:
        pathGTDriftGlobalResults + "analyses_summaries/BLASTP_results/blastp_results.csv",
        pathGTDriftGlobalResults + "sorted_taxonomy.csv"
    output:
        pathGTDriftGlobalResults + "analyses_summaries/table_results/krab_data.csv",
        pathGTDriftGlobalResults + "analyses_summaries/table_results/krabzf_data.csv",
        pathGTDriftGlobalResults + "analyses_summaries/table_results/zf_count.csv"
    shell:
        "python3 "+pathGTDriftScripts+"/analyses/prdm9_protein_analysis/python/krab.py -i "+ pathGTDriftData + " -o " + pathGTDriftGlobalResults +"\
        && python3 "+pathGTDriftScripts+"/analyses/prdm9_protein_analysis/python/krabzf.py  -i "+ pathGTDriftData + " -o " + pathGTDriftGlobalResults +"\
        && python3 "+pathGTDriftScripts+"analyses/prdm9_protein_analysis/python/zf_analysis.py  -i "+ pathGTDriftData + " -o " + pathGTDriftGlobalResults

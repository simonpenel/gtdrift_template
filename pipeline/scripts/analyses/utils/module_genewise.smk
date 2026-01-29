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

rule filter_genomic:
   """
    Filter genomic sequences.
    """
    input:
        pathGTDriftData+ "genome_assembly/{accession}/genome_seq/genomic.fna"
    output:    
        temp(pathGTDriftData+ "genome_assembly/{accession}/genome_seq/genomic.no_lc.fna")
    shell:
        """
        python3 ../utils/python/remove_low_compl.py {input} {output}
        """
        
rule get_blast_db_no_lc:
    """
    Generate a parsable blast db with no low complexity.
    """
    input:
        fasta = pathGTDriftData+ "genome_assembly/{accession}/genome_seq/genomic.no_lc.fna"
    output:
        temp(directory("data/blastdb_nucleotide_seq/{accession}/db_no_lc/"))        
    shell:
        """
        #mkdir data/blastdb_nucleotide_seq/{wildcards.accession}/db_no_lc &&
        makeblastdb -max_file_sz 4G -in {input.fasta}  -out data/blastdb_nucleotide_seq/{wildcards.accession}/db_no_lc/nucldb -dbtype nucl -parse_seqids
        """

rule get_blast_db_lc:
    """
    Generate a parsable blast db in wich low complexity is not removed.
    """
    input:
        fasta = pathGTDriftData+ "genome_assembly/{accession}/genome_seq/genomic.fna"
    output:
        temp(directory("data/blastdb_nucleotide_seq/{accession}/db_lc/"))        
    shell:
        """
        #mkdir data/blastdb_nucleotide_seq/{wildcards.accession}/db_lc &&
        makeblastdb -max_file_sz 4G -in {input.fasta}  -out data/blastdb_nucleotide_seq/{wildcards.accession}/db_lc/nucldb -dbtype nucl -parse_seqids
        #formatdb -i {input} -n data/blastdb_nucleotide_seq/{wildcards.accession}/nucldb -p F -o T
        """
        

rule get_tblastn:
    """
    Run tblastn on database on db with no complexity.
    """
    input:
        cds=pathGTDriftResource+"ref_align/Prdm9_Metazoa_Reference_alignment/exon_peptides/{exon}.fst",        
        db="data/blastdb_nucleotide_seq/{accession}/db_no_lc/"
    output:
        "results/{accession}/Step1_blast/tblastn/PRDM9_{exon}.tblastn.fmt7"
    shell:
        """
        /beegfs/home/penel/tmpdir/ncbi/ncbi-blast-2.16.0+-src/c++/ReleaseMT/bin/tblastn   -query {input.cds} -db data/blastdb_nucleotide_seq/{wildcards.accession}/db_no_lc/nucldb -out {output} -evalue 1e-3 -max_target_seqs 500 -max_hsps 180 -outfmt "7 delim=  qseqid qlen sseqid slen pident nident length mismatch gapopen qstart qend sstart send bitscore evalue" -num_threads 4
        """

def exon_done_tblastn(wildcards):
    return expand("results/" + wildcards.accession + "/Step1_blast/tblastn/PRDM9_{exon}.tblastn.fmt7", exon=EXONS)
    
rule get_loci:
    """
    Get all candidate loci.
    """
    input:
        exon_done_tblastn
    output:
        "results/{accession}/Step1_blast/tblastn/candidate_loci"
    shell:
        """
        python3 ../utils/python/get_loci.py -i results/{wildcards.accession}/Step1_blast/tblastn -o {output} -t nucl
        """
        
rule get_chromosomes:
    """
    Prepare a blastdbcmd batch entry file for extraction.
    """
    input:
        "results/{accession}/Step1_blast/tblastn/candidate_loci"
    output:
        "results/{accession}/Step2_extract_loci/separated_candidates.txt"
    shell:
        """
        python3 ../utils/python/get_chromosomes.py -i {input} -o results/{wildcards.accession}/Step2_extract_loci/separated_candidates.txt
        """
        
rule extract_sequences:
    """
    Extract candidate loci.
    """
    input:
        entry="results/{accession}/Step2_extract_loci/separated_candidates.txt",
        db="data/blastdb_nucleotide_seq/{accession}/db_lc/"
    output:
        "results/{accession}/Step2_extract_loci/candidate_loci/candidates.fna"
    shell:
        """
        blastdbcmd -db data/blastdb_nucleotide_seq/{wildcards.accession}/db_lc/nucldb -entry_batch {input.entry} > {output}
        """
       
rule annotate_candidates:
    """
    Giving each candidate a unique ID for later parsing
    """
    input:
        "results/{accession}/Step2_extract_loci/candidate_loci/candidates.fna"
    output:
        "results/{accession}/Step2_extract_loci/candidate_loci/annotated_candidates.fna"
    shell:
        """
        python3 ../utils/python/annotate_candidates.py -i {input} -o {output}
        """  
checkpoint separate_sequences:
    """
    Separate candidate loci for genewise analysis.
    """
    input:
        "results/{accession}/Step2_extract_loci/candidate_loci/annotated_candidates.fna"
    output:
        directory("results/{accession}/Step2_extract_loci/separated_candidates")
    shell:
        """
        mkdir {output}\
        && python3 ../utils/python/separate_candidates.py -i {input} -o {output}
        """

def aggregate(wildcards):
    candidate_output = checkpoints.separate_sequences.get(accession=wildcards.accession).output[0]
    return expand("results/{accession}/Step2_extract_loci/separated_candidates/{candidate}.fna", accession=wildcards.accession, candidate=glob_wildcards(os.path.join(candidate_output, f"{{candidate}}.fna")).candidate)

rule genewisedb:
    """
    Genewisedb on candidate loci.
    Snakemake is very tedious when it comes to jobs with an unknown number of outputs.
    Therefore the shell script here is rather long. Perhaps it can be made into a separate file?
    """
    input:
        candidates=aggregate,
        ref=pathGTDriftResource+"ref_align/Prdm9_Metazoa_Reference_alignment/PRDM9_metazoa_ReferenceSequences.fa"

    output:
        touch(".processed/genewisedb.{accession}.done"),
        res=".processed/genewisedb.{accession}.lst"
    shell:
        """
        accession={wildcards.accession};
        echo ${{accession}};
        mkdir -p results/"${{accession}}"/Step3_genewise ;
        elt=' ' read -a array <<< "{input.candidates}";
        i=1;
        echo "processing genewise on {input.candidates}..."
        for cand in ${{array[@]}};
        do
            file_output=results/"${{accession}}"/Step3_genewise/"${{i}}".gw;
            echo "processing genewise..."
            echo ${{file_output}};
            genewisedb {input.ref} ${{cand}} -prodb -dnas -genes -pseudo -cdna -pep -quiet -init local -subs 1e-6 -indel 1e-6 -pretty -aln {ALN} > ${{file_output}} || touch ${{file_output}};
            python3 ../utils/python/check_end_file.py -i ${{file_output}}
            echo "completed"
            echo ${{file_output}} >> {output.res}
            ((i++))
        done
        touch {output.res}
        """         
rule concatenate_candidates:
    """
    Concatenating candidates for genewise parsing
    """
    input:
        ".processed/genewisedb.{accession}.done",
        res=".processed/genewisedb.{accession}.lst"
    output:
        "results/{accession}/Step3_genewise/gw.concat"
    shell:
        """
        export testgw=`wc -l {input.res}|cut -f1 -d" "`;
        touch {output} 
        echo $testgw
        if [ $testgw != "0" ];
        then
        echo "Some results for genewise"
        cd results/{wildcards.accession}/Step3_genewise;
        cat *.gw > gw.concat;
        cd ../../..;
        #rm genewisedb.{wildcards.accession}.done
        else 
        echo "No results for genewise"
        fi
        """          

rule genewise_parser:
    """
    Parse through generated genewise files
    """
    input:
        concat="results/{accession}/Step3_genewise/gw.concat",
        res=".processed/genewisedb.{accession}.lst"
    output:
        fasta="results/{accession}/Step3_genewise/genewise_predicted_proteins.faa",
        text="results/{accession}/Step3_genewise/genewise_prediction.txt"
    shell:
        """
        export testgw=`wc -l {input.res}|cut -f1 -d" "`;
        touch {output.text}
        touch {output.fasta} 
        echo $testgw
        if [ $testgw != "0" ];
        then
        echo "Some results for genewise"
        python3 ../utils/python/genewise_parser.py -i {input.concat} -a {wildcards.accession} -o {output.text} -f {output.fasta}
        else 
        echo "No results for genewise"
        fi
        """
rule hmm_search:
    """
    Proteome search using the HMMs.
    """
    input:
        model=get_reference_file,
        protein="results/{accession}/Step3_genewise/genewise_predicted_proteins.faa",
    output:
        # table of of per-sequence hits
        table=pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/tbl/{domain}",
        # table of per-domain hits
        domains=pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/domtbl/{domain}_domains",
    shell:
        """
        export testgw=`wc -l {input.protein}|cut -f1 -d" "`;
        touch {output.table}
        touch {output.domains} 
        echo $testgw
        if [ $testgw != "0" ];
        then
        {RUNCMD} hmmsearch -E 1E-3 --domE 1E-3 --tblout {output.table} --domtblout {output.domains} --noali {input.model} {input.protein}
        else 
        echo "No candidate"
        fi
        """

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
        "../utils/python/hmmsearch_parser.py"

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
        # per-domain hits in tabular format in which overlapping zinc
        # finger domains are merged to create one big domain with 
        # multiple repetitions.
        domain_summary=pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/domtbl/{domain}_domains_summary",
    script:
        "../utils/python/domain_parser.py"
 
# ------------------------------------------------------------
# Function sending the accession number
# ------------------------------------------------------------
def accession_nb(wildcards):
    return wildcards.accession

# ------------------------------------------------------------
# summarize_hmm_results
# write results of hmm search  for each assembly and domain.
# Output format example:
# ;SeqID;SET Query;SET E-value;SET Score;Nb SET domains;SET domain start;SET domain end;Taxid
# ------------------------------------------------------------
rule summarize_hmm_results:
    """
    Creation of a summary table of hmm_search results.
    """
    input:
        # organisms_data file
        organisms_file=pathGTDriftData + "organisms_data",
        # path of all per-sequence hits in tabular format 
        domain_per_sequence_tabulated=pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/tbl/{domain}_tabulated",
        # path of all per-domain hits in tabular format with overlapping zinc finger domains                     
        domain_per_domain_summary=pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/domtbl/{domain}_domains_summary",        
    output:
        # domain protein statistics for each assembly.
        pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +"summary_hmmsearch_{accession}_{domain}.csv",
    params:    
        accession=accession_nb,  
    script:
        "../utils/python/table_domain_builder_single_genewise.py"


checkpoint check_for_error:
    """
    Check for error during genewisedb
    """
    input:
 #       concat_or_not
 #       "results/{accession}/Step2_extract_loci/annotated_predicted_proteins.faa"
        "results/{accession}/Step3_genewise/genewise_predicted_proteins.faa"
    output:
        directory("results/{accession}/error_check")
    shell:
        """
        mkdir {output}\
        && python3 ../utils/python/check_for_error.py -i {input} -o {output}
        """

def error_check(wildcards):
    err_out = checkpoints.check_for_error.get(**wildcards).output[0]
    ERR = glob_wildcards(os.path.join(err_out, "errcheck.{err}")).err[0]
    return 'results/' + wildcards.accession + '/error_check/errcheck.' + ERR
            
def domain_done(wildcards):
    return expand(pathGTDriftData
        + "genome_assembly/" + wildcards.accession + "/analyses/" + GENOME_RESULTS + "hmm_search/domtbl/{domain}_domains_summary" , domain=ALL_DOMAINS)
        

rule table_editing:
    """
    Creation of a summary table of hmm_search results.
    """
    input:
        domain_done,
        err = error_check
    output:
        temp("results/{accession}/Result_tables/summary_table_prdm9_{accession}.csv")
        #"results/{accession}/Result_tables/summary_table_prdm9_{accession}.csv"      
    shell:
        "python3 ../utils/python/table_builder_genewise.py -a {wildcards.accession} -e {input.err} -o {output}"

rule candidate_parsing:
    """
    Parsing candidates based on several criteria, keeping one per locus.
    """
    input:
        pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +"whole_summary_genewise.csv",
        #"results/{accession}/Result_tables/summary_table_prdm9_{accession}.csv"
    output:
       pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +"parsed_whole_summary_genewise.csv",
    shell:
        """
        echo python3 python/candidate_parser_genewise.py -i {input} -o {output}
        python3 ../utils/python/candidate_parser_genewise.py -i {input} -o {output}
        """

rule select_candidates:       
    input:
        pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +"parsed_whole_summary_genewise.csv",
    output:
        pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +"wad_parsed_whole_summary_genewise.csv"
    shell:
        """
        echo python3 python/candidate_parser_genewise.py -i {input} -o {output}
        python3 ../utils/python/candidate_selection_genewise.py -i {input} -o {output}
        """     

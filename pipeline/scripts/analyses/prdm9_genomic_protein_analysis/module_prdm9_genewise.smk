configfile: "config.json"

if config["mode"] == "guix":
    RUNCMD="guix shell hmmer -- "
else:
    RUNCMD=""
    
storagetype  = config["storagetype"]
#MINI = config["mini"]


if config["Nb_aln_genewise"] == "":
    ALN = '10'
else:
    ALN = config["Nb_aln_genewise"]

# Function to load JSON files
def load_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

# Assign environment variables
globals().update(load_json("../environment_path.json"))


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
#            iget -f /lbbeZone/home/penel/gtdrift/genome_seq/$genomic {pathGTDriftData}"genome_assembly/{wildcards.accession}/genome_seq/$genomic"
#            echo ln -s {pathGTDriftData}"genome_assembly/{wildcards.accession}/genome_seq/$genomic {output.fasta}"
#            ln -s {pathGTDriftData}genome_assembly/{wildcards.accession}/genome_seq/$genomic {output.fasta}
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
        python3 python/remove_low_compl.py {input} {output}
        """

rule get_blast_db:
    """
    Generate a parsable blast db for uncurated sequences.
    """
    input:
 #       fasta="data/assemblies/{accession}/genomic.fna"
        fasta = pathGTDriftData+ "genome_assembly/{accession}/genome_seq/genomic.no_lc.fna",
 #       filter = pathGTDriftData+ "genome_assembly/{accession}/genome_seq/genomic.ans"
    output:
        #temp("data/blastdb_nucleotide_seq/{accession}/nucldb.nhr"),
        #temp("data/blastdb_nucleotide_seq/{accession}/nucldb.nin"),
        #"data/blastdb_nucleotide_seq/{accession}/nucldb.nog"),
        #temp("data/blastdb_nucleotide_seq/{accession}/nucldb.nsd"),
        #temp("data/blastdb_nucleotide_seq/{accession}/nucldb.nsi"),
        #temp("data/blastdb_nucleotide_seq/{accession}/nucldb.nsq")
        temp(directory("data/blastdb_nucleotide_seq/{accession}/db/"))        
    shell:
        """
        mkdir data/blastdb_nucleotide_seq/{wildcards.accession}/db
        makeblastdb -max_file_sz 4G -in {input.fasta}  -out data/blastdb_nucleotide_seq/{wildcards.accession}/db/nucldb -dbtype nucl -parse_seqids
        #formatdb -i {input} -n data/blastdb_nucleotide_seq/{wildcards.accession}/nucldb -p F -o T
        """

rule get_blast_db_curated:
    input:
        #fasta = expand("data/assemblies/{accession}/protein.faa", accession = CURATED)
        fasta = expand(pathGTDriftData+ "genome_assembly/{accession}/annotation/protein.faa", accession = CURATED)        
    output:
        expand("data/blastdb_protein_seq/{accession}/protdb.phr", accession = CURATED),
        expand("data/blastdb_protein_seq/{accession}/protdb.pin", accession = CURATED),
        #expand("data/blastdb_protein_seq/{accession}/protdb.pog", accession = CURATED),
        expand("data/blastdb_protein_seq/{accession}/protdb.psd", accession = CURATED),
        expand("data/blastdb_protein_seq/{accession}/protdb.psi", accession = CURATED),
        expand("data/blastdb_protein_seq/{accession}/protdb.psq", accession = CURATED)
    shell:
        """
        elt=' ' read -r -a array <<< "{CURATED}";
        accession="{CURATED}";
        echo $array;
        for fasta in ${{array[@]}};
        do
            mkdir -p data/blastdb_protein_seq/"$fasta";
            makeblastdb  -max_file_sz 3G  -in data/assemblies/"$fasta"/protein.faa -out data/blastdb_protein_seq/"$fasta"/protdb -dbtype prot -parse_seqids;
            #formatdb -i data/assemblies/"$fasta"/protein.faa -n data/blastdb_protein_seq/"$fasta"/protdb -p T -o T;

        done
        """

rule get_tblastn:
    """
    Run tblastn on database using refs.
    """
    input:
        #cds="data/ref_align/Prdm9_Metazoa_Reference_alignment/exon_peptides/{exon}.fst",
        cds=pathGTDriftResource+"ref_align/Prdm9_Metazoa_Reference_alignment/exon_peptides/{exon}.fst",        
        #nhr="data/blastdb_nucleotide_seq/{accession}/nucldb.nhr",
        #nin="data/blastdb_nucleotide_seq/{accession}/nucldb.nin",
        ##nog="data/blastdb_nucleotide_seq/{accession}/nucldb.nog",
        #nsd="data/blastdb_nucleotide_seq/{accession}/nucldb.nsd",
        #nsi="data/blastdb_nucleotide_seq/{accession}/nucldb.nsi",
        #nsq="data/blastdb_nucleotide_seq/{accession}/nucldb.nsq"
        db="data/blastdb_nucleotide_seq/{accession}/db/"
    output:
        "results/{accession}/Step1_blast/tblastn/PRDM9_{exon}.tblastn.fmt7"
    shell:
        """
 #       tblastn  -db_soft_mask 11 -query {input.cds} -db data/blastdb_nucleotide_seq/{wildcards.accession}/nucldb -out {output} -evalue 1e-3 -max_target_seqs 500 -max_hsps 180 -outfmt "7 delim=  qseqid qlen sseqid slen pident nident length mismatch gapopen qstart qend sstart send bitscore evalue" -num_threads 4
        /beegfs/home/penel/tmpdir/ncbi/ncbi-blast-2.16.0+-src/c++/ReleaseMT/bin/tblastn   -query {input.cds} -db data/blastdb_nucleotide_seq/{wildcards.accession}/db/nucldb -out {output} -evalue 1e-3 -max_target_seqs 500 -max_hsps 180 -outfmt "7 delim=  qseqid qlen sseqid slen pident nident length mismatch gapopen qstart qend sstart send bitscore evalue" -num_threads 4
        """

rule get_blastp:
    """
    Run blastp on database using refs if protein is annotated.
    """
    input:
        expand("data/blastdb_protein_seq/{accession}/protdb.phr", accession = CURATED),
        expand("data/blastdb_protein_seq/{accession}/protdb.pin", accession = CURATED),
        #expand("data/blastdb_protein_seq/{accession}/protdb.pog", accession = CURATED),
        expand("data/blastdb_protein_seq/{accession}/protdb.psd", accession = CURATED),
        expand("data/blastdb_protein_seq/{accession}/protdb.psi", accession = CURATED),
        expand("data/blastdb_protein_seq/{accession}/protdb.psq", accession = CURATED),
        cds=expand("data/ref_align/Prdm9_Metazoa_Reference_alignment/exon_peptides/{exon}.fst", exon = EXON)
    output:
        expand("results/{accession}/Step1_blast/blastp/PRDM9_{exon}.blastp.fmt7", accession = CURATED, exon = EXON)
    shell:
        """
        elt=' ' read -r -a array <<< "{CURATED}"
        for prot in ${{array[@]}};
        do
            mkdir -p results/"$prot"/blastp;
            elt=' ' read -r -a cds <<< "{EXON}";
            for exon in ${{cds[@]}};
            do
#                blastp -query data/ref_align/Prdm9_Metazoa_Reference_alignment/exon_peptides/"$exon".fst -db data/blastdb_protein_seq/"$prot"/protdb -out results/"$prot"/Step1_blast/blastp/PRDM9_"$exon".blastp.fmt7 -evalue 1e-3 -max_target_seqs 500 -max_hsps 180 -outfmt "7 delim=  qseqid qlen sseqid slen pident nident length mismatch gapopen qstart qend sstart send bitscore evalue" -num_threads 4
                blastp -query "+pathGTDriftResource+"/ref_align/Prdm9_Metazoa_Reference_alignment/exon_peptides/"$exon".fst -db data/blastdb_protein_seq/"$prot"/protdb -out results/"$prot"/Step1_blast/blastp/PRDM9_"$exon".blastp.fmt7 -evalue 1e-3 -max_target_seqs 500 -max_hsps 180 -outfmt "7 delim=  qseqid qlen sseqid slen pident nident length mismatch gapopen qstart qend sstart send bitscore evalue" -num_threads 4
            done;
        done
        """

def exon_done_tblastn(wildcards):
    return expand("results/" + wildcards.accession + "/Step1_blast/tblastn/PRDM9_{exon}.tblastn.fmt7", exon=EXON)

def exon_done_blastp(wildcards):
    return expand("results/{accession}/Step1_blast/blastp/PRDM9_{exon}.blastp.fmt7", accession=CURATED, exon=EXON)

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
        python3 python/get_loci.py -i results/{wildcards.accession}/Step1_blast/tblastn -o {output} -t nucl
        """

rule get_loci_curated:
    """
    Get all candidate loci for curated prot database.
    """
    input:
        exon_done_blastp,
        #gff = expand("data/assemblies/{accession}/genomic.gff", accession = CURATED)
        gff = expand(pathGTDriftData+ "genome_assembly/{accession}/annotation/genomic.gff", accession = CURATED)

    output:
        expand("results/{accession}/Step1_blast/blastp/candidate_loci", accession = CURATED)
    shell:
        """
        elt=' ' read -r -a array <<< "{CURATED}"
        for acc in ${{array[@]}};
        do
            python3 python/get_loci.py -i results/"$acc"/Step1_blast/blastp -o results/"$acc"/Step1_blast/blastp/candidate_loci -t prot -g data/assemblies/"$acc"/genomic.gff
        done
        """

rule extract_loci_curated:
    input:
        expand("data/blastdb_protein_seq/{accession}/protdb.phr", accession = CURATED),
        expand("data/blastdb_protein_seq/{accession}/protdb.pin", accession = CURATED),
        #expand("data/blastdb_protein_seq/{accession}/protdb.pog", accession = CURATED),
        expand("data/blastdb_protein_seq/{accession}/protdb.psd", accession = CURATED),
        expand("data/blastdb_protein_seq/{accession}/protdb.psi", accession = CURATED),
        expand("data/blastdb_protein_seq/{accession}/protdb.psq", accession = CURATED),
        expand("results/{accession}/Step1_blast/blastp/candidate_loci", accession = CURATED)
    output:
        expand("results/{accession}/Step1_blast/blastp/extracted_loci.faa", accession = CURATED)
    shell:
        """
        elt=' ' read -r -a array <<< "{CURATED}"
        for acc in ${{array[@]}};
        do
            python3 python/extract_loci.py -i results/"$acc"/Step1_blast/blastp/candidate_loci -db data/blastdb_protein_seq/"$acc"/protdb -o results/"$acc"/Step1_blast/blastp/extracted_loci.faa
        done
        """

rule convert_to_gw:
    input:
        seq = expand("results/{accession}/Step1_blast/blastp/extracted_loci.faa", accession = CURATED),
        info = expand("results/{accession}/Step1_blast/blastp/candidate_loci", accession = CURATED)
    output:
        out = expand("results/{accession}/Step3_genewise/blastp_prediction_output.txt", accession = CURATED),
        new = expand("results/{accession}/Step2_extract_loci/annotated_candidates.faa", accession = CURATED)
    shell:
        """
        elt=' ' read -r -a array <<< "{CURATED}"
        for acc in ${{array[@]}};
        do
            python3 python/convert_to_genewise.py -i results/"$acc"/Step1_blast/blastp/candidate_loci -s results/"$acc"/Step1_blast/blastp/extracted_loci.faa -n results/"$acc"/Step2_extract_loci/annotated_candidates.faa -a "$acc" -o results/"$acc"/Step3_genewise/blastp_prediction_output.txt
        done
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
        python3 python/get_chromosomes.py -i {input} -o results/{wildcards.accession}/Step2_extract_loci/separated_candidates.txt
        """


rule extract_sequences:
    """
    Extract candidate loci.
    """
    input:
        entry="results/{accession}/Step2_extract_loci/separated_candidates.txt",
        #nhr="data/blastdb_nucleotide_seq/{accession}/nucldb.nhr",
        #nin="data/blastdb_nucleotide_seq/{accession}/nucldb.nin",
        ##nog="data/blastdb_nucleotide_seq/{accession}/nucldb.nog",
        #nsd="data/blastdb_nucleotide_seq/{accession}/nucldb.nsd",
        #nsi="data/blastdb_nucleotide_seq/{accession}/nucldb.nsi",
        #nsq="data/blastdb_nucleotide_seq/{accession}/nucldb.nsq"
        db="data/blastdb_nucleotide_seq/{accession}/db/"
    output:
        "results/{accession}/Step2_extract_loci/candidate_loci/candidates.fna"
    shell:
        """
        blastdbcmd -db data/blastdb_nucleotide_seq/{wildcards.accession}/db/nucldb -entry_batch {input.entry} > {output}
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
        python3 python/annotate_candidates.py -i {input} -o {output}
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
        && python3 python/separate_candidates.py -i {input} -o {output}
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
 #       ref="data/ref_align/Prdm9_Metazoa_Reference_alignment/PRDM9_metazoa_ReferenceSequences.fa"
        ref=pathGTDriftResource+"ref_align/Prdm9_Metazoa_Reference_alignment/PRDM9_metazoa_ReferenceSequences.fa"

    output:
        touch("genewisedb.{accession}.done"),
        res="genewisedb.{accession}.lst"
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
            python3 python/check_end_file.py -i ${{file_output}}
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
        "genewisedb.{accession}.done",
        res="genewisedb.{accession}.lst"
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
        res="genewisedb.{accession}.lst"
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
        python3 python/genewise_parser.py -i {input.concat} -a {wildcards.accession} -o {output.text} -f {output.fasta}
        else 
        echo "No results for genewise"
        fi
        """

rule concatenate_predictions:
    input:
        genew = "results/{accession}/Step3_genewise/genewise_prediction.txt",
        blastp = "results/{accession}/Step3_genewise/blastp_prediction_output.txt"
    output:
        "results/{accession}/Step3_genewise/concatenated_genewise_prediction.txt"
    shell:
        """
        cat {input.genew} {input.blastp} > {output}
        """

def curated_check(wildcards):
    if {wildcards.accession} in CURATED:
        return expand("results/{accession}/Step2_extract_loci/annotated_candidates.faa", accession = CURATED)
    else:
        return expand("results/{accession}/Step3_genewise/genewise_predicted_proteins.faa", accession = ACCESSNB)

rule concatenate_proteins:
    input:
        blastp = curated_check,
        genew = "results/{accession}/Step3_genewise/genewise_predicted_proteins.faa"
    output:
        "results/{accession}/Step3_genewise/annotated_predicted_proteins.faa"
    shell:
        """
        elt=' ' read -r -a array <<< "{input.blastp}"
        for path in ${{array[@]}};
        do
            if [ "$path" == "results/{wildcards.accession}/Step2_extract_loci/annotated_candidates.faa" ];
            then
                cat {input.genew} results/{wildcards.accession}/Step2_extract_loci/annotated_candidates.faa > {output}
            elif [ "$path" == "results/{wildcards.accession}/Step3_genewise/genewise_predicted_proteins.faa" ];
            then
                cp {input.genew} {output}
            fi
        done
        """

def concat_or_not(wildcards):
    if {wildcards.accession} in CURATED:
        return "results/{accession}/Step2_extract_loci/annotated_predicted_proteins.faa"
    else:
        return "results/{accession}/Step3_genewise/annotated_predicted_proteins.faa"

checkpoint check_for_error:
    """
    Check for error during genewisedb
    """
    input:
        concat_or_not
    output:
        directory("results/{accession}/error_check")
    shell:
        """
        mkdir {output}\
        && python3 python/check_for_error.py -i {input} -o {output}
        """

def error_check(wildcards):
    err_out = checkpoints.check_for_error.get(**wildcards).output[0]
    ERR = glob_wildcards(os.path.join(err_out, "errcheck.{err}")).err[0]
    return 'results/' + wildcards.accession + '/error_check/errcheck.' + ERR

rule get_protdb:
    """
    Generate a blast db for the predicted proteins.
    """
    input:
        concat=concat_or_not,
        res="genewisedb.{accession}.lst"
    output:
        phr="results/{accession}/Step3_genewise/protdb.phr",
        pin="results/{accession}/Step3_genewise/protdb.pin",
        #pog="results/{accession}/Step3_genewise/protdb.pog",
        psq="results/{accession}/Step3_genewise/protdb.psq",
        psi="results/{accession}/Step3_genewise/protdb.psi",
        psd="results/{accession}/Step3_genewise/protdb.psd"
    shell:
        """
        export testgw=`wc -l {input.res}|cut -f1 -d" "`;
        touch {output}
        echo $testgw
        if [ $testgw != "0" ];
        then
        makeblastdb  -max_file_sz 3G  -in {input.concat} -title protdb -out results/{wildcards.accession}/Step3_genewise/protdb -dbtype prot -parse_seqids
        #formatdb -i {input} -t protdb -n results/{wildcards.accession}/Step3_genewise/protdb -p T -o T
        else 
        echo "No results for genewise"
        fi
        """
         
rule hmm_build:
    """
    HMM creation.
    """
    input:
 #       "data/ref_align/Prdm9_Metazoa_Reference_alignment/Domain_{domain}_ReferenceAlignment.fa"
        pathGTDriftResource+"ref_align/Prdm9_Metazoa_Reference_alignment/Domain_{domain}_ReferenceAlignment.fa"
    output:
        "data/hmm_build/{domain}.hmm"
    shell:
        "{RUNCMD} hmmbuild {output} {input}"
 
rule hmm_search:
    """
    Proteome search using the HMMs.
    """
    input:
        model="data/hmm_build/{domain}.hmm",
        protein="results/{accession}/Step3_genewise/genewise_predicted_proteins.faa",
        res="genewisedb.{accession}.lst"
    output:
        table = "results/{accession}/Step4_Hmm/tbl/{domain}",
        domains = "results/{accession}/Step4_Hmm/domtbl/{domain}_domains"
    shell:
        """
        export testgw=`wc -l {input.res}|cut -f1 -d" "`;
        touch {output.table}
        touch {output.domains} 
        echo $testgw
        if [ $testgw != "0" ];
        then
        {RUNCMD} hmmsearch -E 1E-3 --domE 1E-3 --tblout {output.table} --domtblout {output.domains} --noali {input.model} {input.protein}
        else 
        echo "No results for genewise"
        fi
        """

rule tbl_processing:
    """
    Result file processing for a later use.
    """
    input:
        "results/{accession}/Step4_Hmm/tbl/{domain}"
    output:
        "results/{accession}/Step4_Hmm/tbl/{domain}_processed"
    shell:
        "python3 python/hmmsearch_parser.py -i {input} -o {output}"

rule domain_processing:
    """
    Result file processing for a later use.
    """
    input:
        "results/{accession}/Step4_Hmm/tbl/{domain}_processed",
        domain_data="results/{accession}/Step4_Hmm/domtbl/{domain}_domains"
    output:
        processed="results/{accession}/Step4_Hmm/domtbl/{domain}_domains_processed",
        summary="results/{accession}/Step4_Hmm/domtbl/{domain}_domains_summary"
    shell:
       "python3 python/domain_parser.py -i {input.domain_data} -o {output.processed} -s {output.summary}"


def domain_done(wildcards):
    return expand("results/" + wildcards.accession + "/Step4_Hmm/domtbl/{domain}_domains_summary" , domain=DOMAIN)

rule table_editing:
    """
    Creation of a summary table of hmm_search results.
    """
    input:
        domain_done,
        err = error_check
    output:
        temp("results/{accession}/Result_tables/summary_table_prdm9_{accession}.csv")
    shell:
        "python3 python/table_builder_genewise.py -a {wildcards.accession} -e {input.err} -o {output}"

rule candidate_parsing:
    """
    Parsing candidates based on several criteria, keeping one per locus.
    """
    input:
        "results/{accession}/Result_tables/summary_table_prdm9_{accession}.csv"
    output:
        "results/{accession}/Result_tables/parsed_summary_table_prdm9_{accession}.csv"
    shell:
        """
        python3 python/candidate_parser_genewise.py -i {input} -o {output}
        """

rule read_table:
    """
    Reads each summary table and runs a blastp analysis on every candidate
    """
    input:
        table="results/{accession}/Result_tables/parsed_summary_table_prdm9_{accession}.csv",
        phr="results/{accession}/Step3_genewise/protdb.phr",
        pin="results/{accession}/Step3_genewise/protdb.pin",
        #pog="results/{accession}/Step3_genewise/protdb.pog",
        psq="results/{accession}/Step3_genewise/protdb.psq",
        psi="results/{accession}/Step3_genewise/protdb.psi",
        psd="results/{accession}/Step3_genewise/protdb.psd",
        prdmfam=pathGTDriftResource+"PRDM_family_HUMAN"
    output:
        "results/{accession}/Step4_Hmm/blastp.txt",
        "results/{accession}/Result_tables/summary_table_{accession}.csv"
    shell:
        """
        python3 python/blastp_analysis_genewise.py {input.table} {wildcards.accession} {input.prdmfam}\
        """

rule summary:
    """
    Concatenation of each proteome blastp results.
    """
    input: 
        expand("results/{accession}/Step4_Hmm/blastp.txt", accession=ACCESSNB)
    output:
        "results/BLASTP_results/blastp_summary.txt"
    shell:
        """
        cat {input} > {output}
        """

rule blastp_results:
    """
    Writing a table from the concatenation
    """
    input:
        "results/BLASTP_results/blastp_summary.txt"
    output:
        "results/BLASTP_results/blastp_results.csv",
    shell:
        """
        python3 python/blastp_table.py\
        """

rule taxonomy:
    """
    Creation of a table associating a genome accession number to its complete taxonomy
    """
    input:
        "results/BLASTP_results/blastp_summary.txt"
    output:
        "data/resources/sorted_taxonomy.csv"
    shell:
        "python3 python/taxonomy.py"

def domain_processed(wildcards):
    return expand("results/{accession}/Step4_Hmm/tbl/{domain}_processed" , domain=DOMAIN,accession=ACCESSNB)


rule create_table:
    """
    Creation of multiple result table using blastp results and hmm search results
    """
    input:
        "results/BLASTP_results/blastp_results.csv",
        "data/resources/sorted_taxonomy.csv"
    #    domain_processed
    output:
        pathGTDriftGlobalResults+"prdm9_genomic_protein_analysis/summarized_results/krab_data.csv",
        pathGTDriftGlobalResults+"prdm9_genomic_protein_analysis/summarized_results/krabzf_data.csv",
        pathGTDriftGlobalResults+"prdm9_genomic_protein_analysis/summarized_results/zf_count.csv",
        pathGTDriftGlobalResults+"prdm9_genomic_protein_analysis/summarized_results/table_prdm9.csv",

    shell:
        """
        python3 python/krab.py\
        && python3 python/krabzf.py\
        && python3 python/zf_analysis.py\
        && python3 python/table_prdm9.py
        """

rule result_database:
    """
    Creation of a database for the resulting proteins.
    """
    input:
        rules.genewise_parser.output[0]
    output:
        phr="results/{accession}/output_db/out_prot.phr",
        pin="results/{accession}/output_db/out_prot.pin",
        #pog="results/{accession}/output_db/out_prot.pog",
        psq="results/{accession}/output_db/out_prot.psq",
        psi="results/{accession}/output_db/out_prot.psi",
        psd="results/{accession}/output_db/out_prot.psd"
    shell:
        """
        makeblastdb  -max_file_sz 3G  -in {input} -title out_prot -out results/{wildcards.accession}/output_db/out_prot -dbtype prot -parse_seqids
        #formatdb -i {input} -t out_prot -n results/{wildcards.accession}/output_db/out_prot -p T -o T 
        """

rule get_result_batch:
    """
    Entry batch of relevant hits
    """
    input:
        "results/{accession}/summary_table_{accession}.csv"
    output:
        "results/{accession}/predict_batch.txt"
    shell:
        """
        python3 python/get_predictions.py -i {input} -o {output}
        """

rule extract_proteins:
    """
    Get relevant proteins
    """
    input:
        batch="results/{accession}/predict_batch.txt",
        phr="results/blastdb_protein_seq/{accession}/protdb.phr",
        pin="results/blastdb_protein_seq/{accession}/protdb.pin",
        #pog="results/blastdb_protein_seq/{accession}/protdb.pog",
        psq="results/blastdb_protein_seq/{accession}/protdb.psq",
        psi="results/blastdb_protein_seq/{accession}/protdb.psi",
        psd="results/blastdb_protein_seq/{accession}/protdb.psd"
    output:
        hits="results/{accession}/predicted_proteins/homologue_hits.faa"
    shell:
        """
        blastdbcmd -db results/blastdb_protein_seq/{wildcards.accession}/protdb -entry_batch {input.batch} > {output.hits}
        """

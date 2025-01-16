import json
import os
import pandas as pd

rule build_hmm_for_paralogs:
    input:
        alignments_dir = pathGTDriftResource + "ref_align_for_paralogy_check/" + DOMAIN_HMM_DIR +"{domain}",
    output:
        pathGTDriftResource + "hmm_build/paralogy_check/" + DOMAIN_HMM_DIR +"profil_{domain}.hmm.h3f",
        pathGTDriftResource + "hmm_build/paralogy_check/" + DOMAIN_HMM_DIR +"profil_{domain}.hmm.h3i",
        pathGTDriftResource + "hmm_build/paralogy_check/" + DOMAIN_HMM_DIR +"profil_{domain}.hmm.h3m",
        pathGTDriftResource + "hmm_build/paralogy_check/" + DOMAIN_HMM_DIR +"profil_{domain}.hmm.h3p",
        hmm_db = pathGTDriftResource + "hmm_build/paralogy_check/" + DOMAIN_HMM_DIR +"profil_{domain}.hmm"       
    shell:
        """
        echo "[build_hmm_for_paralogs] Calculating hmm profiles ..." &&
        for file in `ls /beegfs/banque/peneldb/gtdrift_template/pipeline/resources/ref_align_for_paralogy_check/domain_example/{wildcards.domain}/*fst`; do echo "[build_hmm_for_paralogs] Build hmm on $file";output=`basename $file|cut -f1 -d"."`;echo "[build_hmm_for_paralogs] Output is $output.hmm"; {RUNCMD} hmmbuild -o $output.hmmbuild.log  "/beegfs/banque/peneldb/gtdrift_template/pipeline/resources/ref_align_for_paralogy_check/domain_example/{wildcards.domain}/"$output".hmm" $file; done &&
        echo "[build_hmm_for_paralogs] Concatenate hmm profiles ..." &&
        cat /beegfs/banque/peneldb/gtdrift_template/pipeline/resources/ref_align_for_paralogy_check/domain_example/{wildcards.domain}/*.hmm > {output.hmm_db} &&
        echo "[build_hmm_for_paralogs] Compress hmm database" &&
        {RUNCMD} hmmpress  {output.hmm_db}

        """  

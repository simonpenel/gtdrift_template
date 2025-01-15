import json
import os
import pandas as pd

rule build_hmm_for_paralogs:
    input:
        alignments_dir = pathGTDriftResource + "ref_align_for_paralogy_check/" + DOMAIN_HMM_DIR +"{domain}",
    output:
        hmm_db = pathGTDriftResource + "hmm_build/paralogy_check/" + DOMAIN_HMM_DIR +"profil_{domain}.hmm",
    shell:
        """
        for file in `ls /beegfs/banque/peneldb/gtdrift_template/pipeline/resources/ref_align_for_paralogy_check/domain_example/{wildcards.domain}/*fst`; do echo "build hmm on $file";output=`basename $file|cut -f1 -d"."`;echo $output; {RUNCMD} hmmbuild "/beegfs/banque/peneldb/gtdrift_template/pipeline/resources/ref_align_for_paralogy_check/domain_example/{wildcards.domain}/"$output".hmm" $file; done &
        cat /beegfs/banque/peneldb/gtdrift_template/pipeline/resources/ref_align_for_paralogy_check/domain_example/{wildcards.domain}/*.hmm > {output.hmm_db} &
        {RUNCMD} hmmpress  {output.hmm_db}

        """  

# Télécharger les génomes et les annotation à partir du NCBI

Ce pipeline permet de télécharger les génomes et leurs annotations étant donné une liste d'assemblage, dans `/beegfs/banque/gtdrift/data/genome_assembly/XXassemblyXX/genome_seq` et `/beegfs/banque/gtdrift/data/genome_assembly/``XXassemblyXX``/annotation` respectivement. Des liens symboliques seront créés pour faciliter les analyses. Parfois le téléchargement en simultané de plusieurs fichiers bug.


Le fichier de configuration est sous la forme

```
{
  "storagetype": "irods",
  "assembly_list": [
"GCF_029281585.2" , 
"GCF_029289425.2" , 
"GCF_028885625.2" , 
"GCF_028858775.2" , 
"GCF_028885655.2" , 
"GCF_000001405.40" 
]  
}
```
Le champ "storagetype" indique si la séquence du génome (fichier *.fna) doit être stockéé
localement (_local_) ou sur irods (_irods_).


La commande pour lancer ce pipeline :

``` bash
snakemake -s  collecting_annotations.smk --configfile config.json  --cores 1
```


Une commande pour lancer ce pipeline avec un fichier de conf génére par ./prdm9_genomic_protein_analysis/generate_conf_for_collecting_genome_annotation.py 

``` bash
snakemake -s  collecting_annotations.smk --configfile config_prdm9_genomic.json  --cores 1
```

Un autre commande avec singularity

``` bash
snakemake collect_everything --configfile config.json  --use-singularity --singularity-args "--bind /beegfs/:/beegfs/" --cluster "sbatch -J {params.name} -p {params.partition} -N 1 --ntasks={params.ntasks} --mem={params.mem} -t {params.time} -o {params.out} -e {params.err}" --rerun-incomplete --rerun-triggers mtime -j 100 -n --forceall
```

Créer un  dag file :

``` bash
snakemake collect_everything --configfile config.json --forceall --dag | dot -Tpdf > dag-GTDrift.pdf
```

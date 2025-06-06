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

Le champ "assembly_list" donne la liste des assemblages à téléchager. La valeur de ce champ peut être 
générée en utilisant le script décrit dans le pipeline  get_genomes_list.


La commande pour lancer ce pipeline :

``` bash
snakemake -s  collecting_annotations.smk --configfile config.json  --cores 1
```


Une commande pour lancer ce pipeline avec un fichier de conf génére par ./prdm9_genomic_protein_analysis/generate_conf_for_collecting_genome_annotation.py 

``` bash
snakemake -s  collecting_annotations.smk --configfile config_prdm9_genomic.json  --cores 1
```


Créer un  dag file :

``` bash
snakemake --configfile config.json --forceall --dag | dot -Tpdf > dag-GTDrift.pdf
```

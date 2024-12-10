# Lancer la recherche de PRDM9 dans les protéomes

Ce pipeline lance la recherche de PRDM9 dans les protéomes.


Le fichier de configuration est sous la forme

```
{
  "mode": "guix",
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
Le champ "mode" indique les calculs hmm doivent utiliser l'executable local (_local_) ou l'executable via guix (_guix_).

Le champ "assembly_list" donne la liste des assemblages à téléchager. La valeur de ce champ peut être 
générée en utilisant le script décrit dans le pipeline  get_genomes_list.


Lancer l'analyse avec la commande:

```
snakemake -s process_stats_prdm9.smk --cores 1
```

Genere un DAG avec la commande:
```
snakemake -s process_stats_prdm9.smk  --configfile config_dag_single.json --forceall --dag | dot -Tpdf > dag_single_process_stats_prdm9.pdf

snakemake -s process_stats_prdm9.smk  --configfile config_dag_2.json --forceall --dag | dot -Tpdf > dag_2_process_stats_prdm9.pdf

```

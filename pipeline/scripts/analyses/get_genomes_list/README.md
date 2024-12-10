# Recupérerer une liste d'assemblages sur le NCBI

Ce pipeline permet de recupérer la liste des assemblages pour un taxon donné.


La commande pour lancer ce pipeline :

``` bash
snakemake -s  get_list.smk  --cores 1
```

Le  fichier de configuration config.json:
```
{    
    "query": '"\"Hominidae\"[Organism]"'
}
```

Cela va générer le fichier 

```
data/resources/organisms_data 
```
qui peut être copié dans pathGTDriftData/data_results_per_assembly/

# Generating the json files for the prdm9 analysis pipelines

```
python3 generate_json_and_query.py  data/resources/organisms_data  hominidae_assemblies hominidae_queries
```


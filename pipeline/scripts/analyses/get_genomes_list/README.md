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

Dans cet exemple de demonstration, on va récupérer tous les assemblages d'hominidés.

Cela va générer le fichier suivant:

```
data/resources/organisms_data 
```

Lorsque l'on récupère tous les assemblages de métazoaires, ce fichier
 doit être copié dans pathGTDriftData/data_results_per_assembly/, il  servira de référence
 dans différentes étapes des autre pipelines.
  

# Générer les fichier json nécessaires aux autres pipelines


Le fichier de configuration de prdm9_protein_analysis contient une liste d'assemblages versionés au format json.

Le fichier de configuration de prdm9_genomic_protein_analysis contient une requete sur des assemblages non versionés au format json.
 
 Le script generate_json_and_query.py permet de générer ces morceaux de code json.
 

```
python3 generate_json_and_query.py  data/resources/organisms_data  hominidae_assemblies hominidae_queries
```


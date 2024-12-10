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
 
Dans cet exemple avec les hominidés, on lance 
```
python3 generate_json_and_query.py  data/resources/organisms_data  hominidae_assemblies hominidae_queries
```

Ceci génere les 2 fichiers

```
hominidae_queries.col
hominidae_assemblies.col
```

hominidae_assemblies.col contient:
```
  "assembly_list": [
"GCF_029289425.2" , 
"GCF_029281585.2" , 
"GCF_028885625.2" , 
"GCF_028858775.2" , 
"GCF_028885655.2" , 
"GCA_963575185.1" , 
"GCF_000001405.40" 
```

hominidae_queries.col contient:
```
(GCF_029289425) 
 OR (GCF_029281585) 
 OR (GCF_028885625) 
 OR (GCF_028858775) 
 OR (GCF_028885655) 
 OR (GCA_963575185) 
 OR (GCF_000001405) 
```

La redondance entre assemblages GCA et GCF dans organisms_data a été supprimée:

```
 grep panis data/resources/organisms_data 
Pan paniscus	9597	GCA_029289425.3	False	False	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/289/425/GCA_029289425.3_NHGRI_mPanPan1-v2.0_pri/GCA_029289425.3_NHGRI_mPanPan1-v2.0_pri
Pan paniscus	9597	GCF_029289425.2	True	True	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/289/425/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri
```

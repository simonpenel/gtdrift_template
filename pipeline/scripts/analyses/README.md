
# Les étapes à suivre:

## 1. Définir l'environnement dans lequel les calculs sont effectués:

Le fichier _environment_path.json_ est utilisé pour définir l'organisation des répertoires:

Modifier ce fichier en remplaçant _my_directory_ par le répertoire dans lequel se trouve le répertoire _gtdrift_template_.

```json
{
  "pathGTDriftData": "/my_directory/gtdrift_template/data_results_per_assembly/",
  "pathGTDriftGlobalResults": "/my_directory/gtdrift_template/global_results/",
  "pathGTDriftResource": "/my_directory/gtdrift_template/pipeline/resources/",
  "pathGTDriftScripts": "/my_directory/gtdrift_template/pipeline/scripts/",
  "pathGTDriftComputEnv": "/my_directory/gtdrift_template/pipeline/computing_environments/",
  "pathGTDriftLog": "/my_directory/gtdrift_template/temp/log/"
}
```

Les 4 répertoires suivants doivent impérativement être définis:
  * pathGTDriftData : contient 
      * le fichier  _organisms_data_ : description de assemblages
      > Ce fichier est généré par le pipeline snakemake _get_genomes_list/get_list.smk_
      * le répertoire _genome_assembly_ : contient pour chaque assemblage :
          * le répertoire _genome_seq_ : contient le fichier .fna du génome ou son chemin sur iRODS
          > Ce répertoire est rempli par le pipeline snakemake _collecting_genome_annotation/collecting_annotations.smk_
          * le repertoire _annotation_ : contient les fichiers d'annotation *gff, les fichiers fasta des protéines et des cds. S'il le génome n'est pas annoté, les fichiers contiennent un message spécifiant qu'il n'existe pas d'annotation. 
          > Ce répertoire est rempli par le pipeline snakemake _collecting_genome_annotation/collecting_annotations.smk_ 
          * le repertoire _analyse_ : contient le résultat des différentes analyses.
           
  * pathGTDriftGlobalResults : contient les résultats globaux pour chaque analyse.
  * pathGTDriftResource : contient des données utiles pour les analyses, entre autres:
      * le répertoire _ref_align/Prdm9_Metazoa_Reference_alignment/_ qui contient :
          * les alignements de réferences  utilisés pour le calcul des hmm.
          * le répertoire _exon_peptide_ qui contient les fichiers fasta des exons de PRDM9 chez les métazoaires. 
      * le répertoire _hmm_build_ qui contient les profils hmm (calculés à partir de  _ref_align/Prdm9_Metazoa_Reference_alignment/_ lors des analyses)
      * le répertoire _PRDM_family_HUMAN_ qui contient le fichier fasta PRDM_family_HUMAN.fa de la famille PRDM chez les métazoaires (ainsi que la base de données blast  calculée lors des analyses)
  * pathGTDriftScripts : le répertoire des différents pipelines.
  


## 2. Recupérerer la description des assemblages pour un taxon donné

Se déplacer dans le répertoire  _get_genomes_list_ et lancer le pipeline snakemake _get_list.smk_.

La commande pour lancer ce pipeline :

``` bash
snakemake -s  get_list.smk  --cores 1
```

Le taxon est défini dans le fichier de configuration config.json:

```json
{    
    "query": "\"Metazoa\"[Organism]"
}
```

Cela va génerer le fichier  _data/resources/organisms_data_   qui devra être copié dans le répertoire _pathGTDriftData_ pour servir de référence.

> Le fichier _data/resources/organisms_data_ sera copié et utilisé comme référence dans le cas des métazaoires. Mais il est possible de générer un fichier dédié à des jeux de données plus réduits pour des test (sur les hominidés par exemple) auquel cas on ne le copira pas.
  
Le fichier de configuration config.json pour les hominidés:

```json
{    
    "query": "\"Hominidae\"[Organism]"
}
```

## 3. Générer la liste des assemblages au format  fichiers json 

Lancer le  script _python generate_json_and_query.py_ qui se trouve dans le répertoire _get_genomes_list_.

>Cela va générer des fichiers json utiles pour collecter les données et lancer les analyses.

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
```json
  "assembly_list": [
"GCF_029289425.2" , 
"GCF_029281585.2" , 
"GCF_028885625.2" , 
"GCF_028858775.2" , 
"GCF_028885655.2" , 
"GCA_963575185.1" , 
"GCF_000001405.40" 
]
```

hominidae_queries.col contient:
```json
(GCF_029289425) 
 OR (GCF_029281585) 
 OR (GCF_028885625) 
 OR (GCF_028858775) 
 OR (GCF_028885655) 
 OR (GCA_963575185) 
 OR (GCF_000001405) 
```

La redondance entre assemblages GCA et GCF dans organisms_data a été supprimée. Si GCA et GCF sont présent on garde GCF:

```
 grep panis data/resources/organisms_data 
Pan paniscus	9597	GCA_029289425.3	False	False	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/289/425/GCA_029289425.3_NHGRI_mPanPan1-v2.0_pri/GCA_029289425.3_NHGRI_mPanPan1-v2.0_pri
Pan paniscus	9597	GCF_029289425.2	True	True	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/289/425/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri
```

Pour l'analyse _prdm9_protein_analysis_, on doit se limiter aux assemblages pour lesquels il existe une annotation. (On peut aussi vouloir restreindre l'analyse prdm9_genomic_protein_analysis aux assemblages avec annotation)

Pour cela on lance le script avec l'option "curated":

```
python3 generate_json_and_query.py  data/resources/organisms_data  hominidae_assemblies_with_annot hominidae_queries_with_annot curated
```
Le fichier hominidae_assemblies_with_annot.col contient:

```json
"assembly_list": [
"GCF_029281585.2" , 
"GCF_029289425.2" , 
"GCF_028885625.2" , 
"GCF_028858775.2" , 
"GCF_028885655.2" , 
"GCF_000001405.40" 
]
```

Le fichier hominidae_queries_with_annot.col  contient:

```json
(GCF_029281585) 
 OR (GCF_029289425) 
 OR (GCF_028885625) 
 OR (GCF_028858775) 
 OR (GCF_028885655) 
 OR (GCF_000001405)
```

L'assemblage  GCA_963575185.1
```
Gorilla beringei	499232	GCA_963575185.1	False	False	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/963/575/185/GCA_963575185.1_PGDP_GorBer/GCA_963575185.1_PGDP_GorBer
```

n'est plus présent.


## 4. Récupérer les données de séquence et les annotations (si elles existent)
Ce pipeline permet de télécharger les génomes et leurs annotations étant donné une liste d'assemblage. Des liens symboliques seront créés pour faciliter les analyses. Parfois le téléchargement en simultané de plusieurs fichiers bug.

Se déplacer dans le répertoire _collecting_genome_annotation_. Créer le fichier de configuration à l'aide des fichiers générés en *3*, puis lancer le pipeline snakemake _collecting_annotations.smk_.
>Cela va collecter les données de séquences et d'annotation des assemblages. Ces données seront stockées dans les répertoires _genome_seq_ et _annotation_ de chaque assemblage du répertoire _pathGTDriftData/genome_assembly_


Pour les hominidés, le fichier de configuration est sous la forme:

```json
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
Le champ "storagetype" indique si la séquence du génome (fichier *.fna) doit être stockéé localement (_local_) ou sur irods (_irods_).

Le champ "assembly_list" donne la liste des assemblages à téléchager. La valeur de ce champ peut être générée en utilisant le script décrit en *3*.


La commande pour lancer ce pipeline :

``` bash
snakemake -s  collecting_annotations.smk --configfile config.json  --cores 1
```





## 5. Lancer les pipelines d'analyse  des données
Il existe un repertoire par analyse.
Créer le fichier de configuration à l'aide des fichiers générés en *3*, puis 
lancer le pipeline snakemake qui se trouve dans le répertoire de l'analyse.

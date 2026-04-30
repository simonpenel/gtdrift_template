# Settings:

> [!IMPORTANT]
> Il est nécessaire d'installer les commandes _EDirect_:
> https://www.ncbi.nlm.nih.gov/books/NBK179288/  

> [!IMPORTANT]
> Il est nécessaire d'installer le programme _seqkit_ :https://bioinf.shenwei.me/seqkit/download/#installation
> 
> Sinon il est possible d'utiliser  _seqkit_  via _pixi_

> [!IMPORTANT]
> Il est nécessaire d'installer le programme _gffread_ (pour l'étape _collecting_genome_annotation_):https://github.com/gpertea/gffread
> 
> Sinon il est possible d'utiliser  _gffread_  via _pixi_

> [!NOTE]
> Le résultat des exemples de commandes peuvent être différent de ceux donnés dans cette documentation
> car les assemblages proposés par le NCBI évoluent. 


# Les étapes à suivre:

## 1. Définir l'environnement dans lequel les calculs sont effectués:

Le fichier _environment_path.json_ est utilisé pour définir l'organisation des répertoires.
> [!NOTE]
>Le fichier se trouve dans
>`gtdrift_template/pipeline/scripts/analyses`

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
  * pathGTDriftData : ce répetoire  contient 
      * le fichier  _organisms_data_ : description de assemblages
      > Ce fichier est généré par le pipeline snakemake _get_genomes_list/get_list.smk_
      * le répertoire _genome_assembly_ : contient pour chaque assemblage :
          * le répertoire _genome_seq_ : contient le fichier .fna du génome ou son chemin sur iRODS
          > Ce répertoire est rempli par le pipeline snakemake _collecting_genome_annotation/collecting_annotations.smk_
          * le repertoire _annotation_ : contient les fichiers d'annotation *gff, les fichiers fasta des protéines et des cds. Si le génome n'est pas annoté, les fichiers contiennent un message spécifiant qu'il n'existe pas d'annotation. 
          > Ce répertoire est rempli par le pipeline snakemake _collecting_genome_annotation/collecting_annotations.smk_ 
          * le repertoire _analyse_ : contient le résultat des différentes analyses.
           
  * pathGTDriftGlobalResults : contient les résultats globaux pour chaque analyse.
  * pathGTDriftResource : contient différentes données utiles pour les analyses, par exemple pour PRDM9:
      * le répertoire _ref_align/Prdm9_Metazoa_Reference_alignment/_ qui contient :
          * les alignements de réferences  utilisés pour le calcul des hmm.
          * le répertoire _exon_peptide_ qui contient les fichiers fasta des exons de PRDM9 chez les métazoaires. 
      * le répertoire _hmm_build_ qui contient les profils hmm (calculés à partir de  _ref_align/Prdm9_Metazoa_Reference_alignment/_ lors des analyses)
      * le répertoire _PRDM_family_HUMAN_ qui contient le fichier fasta PRDM_family_HUMAN.fa de la famille PRDM chez les métazoaires (ainsi que la base de données blast  calculée lors des analyses)
  * pathGTDriftScripts : le répertoire des différents pipelines.
  


## 2. Récupérerer la description des assemblages pour un taxon donné


Se déplacer dans le répertoire  _gtdrift_template/pipeline/scripts/analyses/get_genomes_list_ et lancer le pipeline snakemake _get_list.smk_.

La commande pour lancer ce pipeline :

``` bash
snakemake -s  get_list.smk  --cores 1
```

La commande pour lancer ce pipeline avec _uv_ (recommandé) :

``` bash
uv run snakemake -s  get_list.smk  --cores 1
```

>Le fichier pyproject.toml en cas d'utilisation de _uv_:

```yaml
[project]
name = "get-genomes-list"
version = "0.1.0"
description = "Get the assemblies associated to a taxon"
readme = "README.md"
requires-python = ">=3.12"
dependencies = [
    "snakemake>=9.19.0",
]
```
> [!TIP]
>Utilisation de _uv_ :https://uv-introduction-da5541.pages.in2p3.fr/1

Le taxon est défini dans le fichier de configuration config.json, ici pour les hominidés:

```json
{    
    "query": "\"Hominidae\"[Organism]"
}
```

Cela va génerer le fichier  _data/resources/organisms_data_   qui devra être copié dans le répertoire _pathGTDriftData_ pour servir de référence dans la suite du pipeline, par exemple pour connaître l'espèce associée a un numero d'assemblage.

> [!IMPORTANT]
> Pour l'analyse des métazoaires le fichier _data/resources/organisms_data_ doit être copié dans _pathGTDriftData_ (i.e. _gtdrift_template/data_results_per_assembly/_) et il sera utilisé comme référence par les autres scripts. Mais il est possible de générer un fichier dédié à des jeux de données plus réduits (pour ensuite génerer des liste d'assemblages destinées a des analyses spécifiques à un taxon par exemple) auquel cas on ne le copiera pas dans _pathGTDriftData_:


Le fichier de configuration config.json pour les métazoaires:

```json
{    
    "query": "\"Metazoa\"[Organism]"
}
```

## 3. Générer la liste des assemblages au format  fichiers json 

Toujours dans le répertoire _gtdrift_template/pipeline/scripts/analyses/get_genomes_list_,
lancer le  script _python generate_json_and_query.py_.

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



La redondance entre assemblages GCA et GCF dans organisms_data a été supprimée. Si GCA et GCF sont présent on garde GCF:

```
 grep panis data/resources/organisms_data 
Pan paniscus	9597	GCA_029289425.3	False	False	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/289/425/GCA_029289425.3_NHGRI_mPanPan1-v2.0_pri/GCA_029289425.3_NHGRI_mPanPan1-v2.0_pri
Pan paniscus	9597	GCF_029289425.2	True	True	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/289/425/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri
```

Dans le cas d'une analyse des protéomes on doit se limiter aux assemblages pour lesquels il existe une annotation. 

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

L'assemblage  GCA_963575185.1
```
Gorilla beringei	499232	GCA_963575185.1	False	False	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/963/575/185/GCA_963575185.1_PGDP_GorBer/GCA_963575185.1_PGDP_GorBer
```

n'est plus présent.


## 4. Récupérer les données de séquence et les annotations (si elles existent)
Ce pipeline permet de télécharger les génomes et leurs annotations étant donné une liste d'assemblage. Des liens symboliques seront créés pour faciliter les analyses. Parfois le téléchargement en simultané de plusieurs fichiers plante.

Se déplacer dans le répertoire _gtdrift_template/pipeline/scripts/analyses/collecting_genome_annotation_. Créer le fichier de configuration _config.json_ à l'aide du fichier _*.col_ générés en *3* , puis lancer le pipeline snakemake _collecting_annotations.smk_.
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


La commande pour lancer ce pipeline avec _uv_ (recommandé) :

``` bash
uv run snakemake -s  collecting_annotations.smk --configfile config.json  --cores 1
```


La commande pour lancer ce pipeline avec _pixi_ (nécessaire sur le cluster pbil) :

``` bash
pixi run snakemake -s  collecting_annotations.smk --configfile config.json  --cores 1
```

Le fichier _pixi.toml_:
```yaml
[workspace]
authors = ["Simon Penel <simon.penel@univ-lyon1.fr>"]
channels = ["conda-forge", "bioconda"]
name = "collecting_genome_annotation_pixi"
platforms = ["linux-64"]
version = "0.1.0"

[tasks]

[dependencies]
snakemake = ">=9.16.3,<10"
gffread = ">=0.12.7,<0.13"
```


## 5. Lancer les pipelines d'analyse  des données

Plusieurs pipelines d'analyse sont disponibles dans le repertoire `pipeline/scripts/analyses` :


- pour une analyse en domaines sur les protéomes :

  - process_stats_domain.smk

    Il s'agit d'une analyse générique  en domaines, ces domaines étant définis dans le fichier de configuration `analyse.json`.

    Pour une analyse PRDM9 ces domaines sont SET, SSXRD, KRAB, ZF.
  - process_zincfinger.smk

    Il s'agit d'une analyse spécifique à PRDM9
  - process_zincfinger_dna.smk

    Il s'agit d'une analyse spécifique à PRDM9
  - process_SET_tyrosines.smk

    Il s'agit d'une analyse spécifique à PRDM9

- pour une analyse en domaines sur les génomes:
  - process_genewise.smk  
    
    Il s'agit d'une analyse générique en domaines, ces domaines étant définis dans le fichier de configuration `analyse.json`.

    Pour une analyse PRDM9 ces domaines sont SET, SSXRD, KRAB, ZF.
  - process_zincfinger_genewise.smk
    
    Il s'agit d'une analyse spécifique à PRDM9
  - process_zincfinger_dna_genewise.smk
    
    Il s'agit d'une analyse spécifique à PRDM9


Créer un répertoire pour chaque analyse,   le répertoire doit contenir un fichier de configuration `analyse.json` associé à l'analyse et  le fichier d'assemblages `assembly.json`  généré en *3*.

### Exemples d'utilisation

  - Le répertoire `verif` pour l'analyse `process_stats_domain.smk`
    https://github.com/simonpenel/gtdrift_template/tree/master/pipeline/scripts/analyses/verif

  - Le répertoire `verif_genewise` pour l'analyse `process_genewise.smk`

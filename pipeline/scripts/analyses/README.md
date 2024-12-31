
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
      * le répertoire _genome_assembly_ : contient pour chaque assemblage :
          * le repertoire _genome_seq_ : contient le fichier .fna du génome ou son chemin sur iRODS
          * le repertoire _annotation_ : contient les fichiers d'annotation *gff, les fichiers fasta des protéines et des cds. S'il le génome n'est pas annoté, les fichiers contiennent un message spécifiant qu'il n'existe pas d'annotation.  
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

Lancer le pipeline snakemake qui se trouve dans le répertoire _get_genomes_list_.
Cela va génerer le fichier  _organisms_data_   qui devra être copié dans le répertoire _pathGTDriftData_ et servira de référence. 
  
## 3. Générer la liste des assemblages au format  fichiers json 

Lancer le  script python qui se trouve dans le répertoire _get_genomes_list_.
Cela va générer des fichier json utiles pour collecter les données et lancer les analyses.

## 4. Récupérer les données de séquence et les annotations (si elles existent)
Créer le fichier de configuration à l'aide des fichiers générés en *3*, puis 
lancer le pipeline snakemake qui se trouve dans le répertoire _collecting_genome_annotation_.
Cela va collecter les données de séquences et d'annotation des assemblages. Ces données seront stockées dans les répertoires _genome_seq_ et _annotation_ de chaque assemblage du répertoire _pathGTDriftData/genome_assembly_.

## 5. Lancer les pipelines d'analyse  des données
Il existe un repertoire par analyse.
Créer le fichier de configuration à l'aide des fichiers générés en *3*, puis 
lancer le pipeline snakemake qui se trouve dans le répertoire de l'analyse.

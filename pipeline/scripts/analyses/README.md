# Définir l'environnement dans lequel les calculs sont effectués:

Le fichier _environment_path.json_ est utilisé pour définir l'organisation des répertoires:

Modifier ce fichier en remplacer _my_directory_ par le répertoire dans lequel se trouve le répertoire _gtdrift_template_.

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

# Les étapes à suivre:


## 1. Recupérerer la description des assemblages pour un taxon donné

Lancer le pipeline snakemake qui se trouve dans le répertoire _get_genomes_list_
  
## 2. Générer la liste des assemblages au format  fichiers json 

Lancer le  script python qui se trouve dans le répertoire _get_genomes_list_.

## 3. Récupérer les données de séquence et les annotations (si elles existent)


## 4. Lancer les pipelines d'analyse  des données

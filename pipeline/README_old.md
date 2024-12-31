Dans ce dossier on retrouve tout le nécessaire pour reproduire les analyses, avec les différents sous-dossiers suivants.

### Computing environments

Contient les programmes, les images singularity, les environnements python *etc.*

-   **`envpy3.11`** C'est un environnement virtuel python (Python 3.11.4), activable avec la ligne de commande:

``` bash
source /beegfs/banque/gtdrift/pipeline/envpy3.11/bin/activate
```

Cette environement et créé avec :

``` bash
/beegfs/banque/gtdrift/pipeline/logiciels/Python-3.11.4/bin/virtualenv --python="/beegfs/banque/gtdrift/pipeline/logiciels/Python-3.11.4/bin/python3.11" "/beegfs/banque/gtdrift/pipeline/envpy3.11/"
```

-   **`logiciels`** Contient une multitude de programmes, qui sont utilisés pour générer les analyses.

-   **`R`** Contient les library de R qui peuvent être utilisés dans les scripts, en ajoutant par exemple :

``` r
.libPaths(c( "/beegfs/banque/gtdrift/pipeline/R/x86_64-pc-linux-gnu-library/4.2" , .libPaths() ))
```

-   **`singularity_sif`** Contient les images singularity et leur recette respective.

### scripts

Contient les Snakefile, pipeline et scripts générant les analyses.

### Resources

Contient les fichiers nécessaire pour des analyses, qui sont communs à toutes les espèces. Par exemple la base de données BUSCO.

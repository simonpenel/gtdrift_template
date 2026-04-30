

# Télécharger les génomes et les annotation à partir du NCBI

Ce pipeline permet de télécharger les génomes et leurs annotations étant donné une liste d'assemblage, dans `/beegfs/banque/gtdrift/data/genome_assembly/XXassemblyXX/genome_seq` et `/beegfs/banque/gtdrift/data/genome_assembly/``XXassemblyXX``/annotation` respectivement. Des liens symboliques seront créés pour faciliter les analyses. Parfois le téléchargement en simultané de plusieurs fichiers bug.

>[!NOTE]
> Plus d'information sur https://github.com/simonpenel/gtdrift_template/tree/master/pipeline/scripts/analyses
>
> ou
>
> `more ../README.md`
Créer un  dag file :

``` bash
snakemake --configfile config.json --forceall --dag | dot -Tpdf > dag-GTDrift.pdf
```

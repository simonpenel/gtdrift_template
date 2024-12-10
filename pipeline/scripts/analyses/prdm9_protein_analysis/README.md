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

Le champ "assembly_list" donne la liste des assemblages à analyser. La valeur de ce champ peut être 
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


# Problèmes d'annotation absente

Il peut arriver que le protéome soit absent (le fichier de protéine contient alors le message NO_PROTEIN_DATA) alors que l'assemblage est annoté.

En général, il s'agit d'un problème lors du téléchargement. Il convient alors de relancer le pipeline _collecting_genome_annotation_.
Mais avant cela il faut naturellement supprimer ces fichiers problématiques ainsi que l'analyse effectuée :

```
export  gc=`grep GC config.json |cut -f2 -d\"`
for acc in $gc
do
grep NO_PROTEIN_DATA ../../../../data_results_per_assembly/genome_assembly/$acc/annotation/protein.faa
if [ $? -gt 0 ]
then
echo "ok"
else
echo "pb"
echo "efface ../../../../data_results_per_assembly/genome_assembly/$acc/analyses/prdm9_prot/"
rm -r ../../../../data_results_per_assembly/genome_assembly/$acc/analyses/prdm9_prot/
echo "efface  ../../../../data_results_per_assembly/genome_assembly/$acc/annotation"
rm -r  ../../../../data_results_per_assembly/genome_assembly/$acc/annotation
fi
done
```

Une fois que pipeline _collecting_genome_annotation_ a tourné et téléchargé les annotations manquantes, on peut relancer le pipeline _process_stats_prdm9.smk_

Dans certains cas (rares) il s'agit d'une erreur dans fichier organism_data due à une information eronnée donné par le NCBI via _esearch_. Il faut alors virer manuellement l'assemblage du fichier de configuration.
On peut aussi ecrire un script pour genere un nouveau fichier de configuration si on est sur que le probleme ne vient pas du téléchargement.

```
head -3 config.json > config.json.clean
export  gc=`grep GC config.json |cut -f2 -d\"`
for acc in $gc
do
grep NO_PROTEIN_DATA ../../../../data_results_per_assembly/genome_assembly/$acc/annotation/protein.faa
if [ $? -gt 0 ]
then
echo "ok"
echo "\"$acc\" ," >> config.json.clean
else
echo "pb"
echo "efface ../../../../data_results_per_assembly/genome_assembly/$acc/analyses/prdm9_prot/"
rm -r ../../../../data_results_per_assembly/genome_assembly/$acc/analyses/prdm9_prot/
fi
done
echo "]}" >> config.json.clean
```

et virer la dernière virgule dans config.json.clean.
On peut remplacer le fichier de configuration par le fichier généré et  relancer le pipeline _process_stats_prdm9.smk_ .

# Graphes

## Graphe pour 1 assemblage

![Image graphe](https://raw.githubusercontent.com/simonpenel/gtdrift_template/refs/heads/master/pipeline/scripts/analyses/prdm9_protein_analysis/dag_single_process_stats_prdm9.png "Graphe pour un assemblage")


## Graphe pour 2 assemblages

![Image graphe](https://raw.githubusercontent.com/simonpenel/gtdrift_template/refs/heads/master/pipeline/scripts/analyses/prdm9_protein_analysis/dag_2_process_stats_prdm9.png "Graphe pour 2 assemblages")



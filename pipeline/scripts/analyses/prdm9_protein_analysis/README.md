Generer la liste de GC en json a partir des donnees collectees
#export  gc=`grep ">"   ../../../../data_results_per_assembly/genome_assembly/*/annotation/protein.faa |cut -f7 -d"/"|sort -u `
#for acc in $gc; do echo "\"$acc\","; done > acc.json


#Vire les fichiers de protein qui contient THERE_IS_NO_PROTEIN_DATA (cas de pb de telechargement)

#head -3 config.json.before_clean > config.json.clean
#export  gc=`grep GC config.json.before_clean |cut -f2 -d\"`
#for acc in $gc
#do
#echo "\"$acc\" ," >> config.json.clean
#grep NO_PROTEIN_DATA ../../../../data_results_per_assembly/genome_assembly/$acc/annotation/#protein.faa
#if [ $? -gt 0 ]
#then
#echo "ok"
#else
#echo "pb"
#rm ../../../../data_results_per_assembly/genome_assembly/$acc/analyses/prdm9_prot/protdb*
#fi
#done
#echo "]}" >> config.json.clean
#et virer la dernire virgule

Lancer l'analyse avec la commande:

```
snakemake -s process_stats_prdm9.smk --cores 1
```

Genere un DAG avec la commande:
```
snakemake -s process_stats_prdm9.smk  --configfile config_dag_single.json --forceall --dag | dot -Tpdf > dag_single_process_stats_prdm9.pdf

snakemake -s process_stats_prdm9.smk  --configfile config_dag_2.json --forceall --dag | dot -Tpdf > dag_2_process_stats_prdm9.pdf

```

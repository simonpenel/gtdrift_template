# Lancer la recherche de de plusieurs domaines dans les protéomes

Ce pipeline lance la recherche de plusieurs dans les protéomes.


Le fichier de configuration est sous la forme

```json

{ 
  "mode": "",
  "analyse_dir_name": "domain_example/",
  "domain_references" : {
"SET":"Domain_SET_ReferenceAlignment2024",
"KRAB":"Domain_KRAB_ReferenceAlignment2024",
"SSXRD":"Domain_SSXRD_ReferenceAlignment2024",
"ZF":"Domain_ZF_unit_ReferenceAlignment2024",
},
  "domains" :  ["SET"],    
  "domains_simple" : ["KRAB","SSXRD","ZF"],   
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

- Le champ "mode" indique si les calculs hmm doivent utiliser l'executable via guix (_guix_).

- Le champ "assembly_list" donne la liste des assemblages à analyser. La valeur de ce champ peut être 
générée en utilisant le script décrit dans le pipeline  get_genomes_list (cf https://github.com/simonpenel/gtdrift_template/tree/master/pipeline/scripts/analyses#3-g%C3%A9n%C3%A9rer-la-liste-des-assemblages-au-format--fichiers-json).


> Lancer l'analyse avec la commande:

```
snakemake -s process_stats_domain.smk -c 1
```



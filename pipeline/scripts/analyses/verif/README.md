PROTEIN DOMAIN ANALYSIS ON PROTEOMES
====================================

## Settings 

We will use "uv" to run snakemake because it is the easisest way of managing python pacakges and versions especially when we hav no root permission.

Type :
```
uv init
uv add pandas
uv add hmmer
uv add Bio
```

R has to be installed with the package Biostrings

## Configuration files

- `analyse.json` file :

```
{ 
  "mode": "",  
  "storagetype": "irods",
  "analyse_dir_name": "verif/",
  "resources_dir_name": "protein_domain_PRDM9/", 
  "domain_references" : {
"SET":"Domain_SET_ReferenceAlignment2024",
"KRAB":"Domain_KRAB_ReferenceAlignment2024",
"SSXRD":"Domain_SSXRD_ReferenceAlignment2024",
"ZF":"Domain_ZF_unit_ReferenceAlignment2024"
},
  "domains" :  ["SET"],    
  "domains_simple" : ["SSXRD","ZF","KRAB"],
  "domain_aln_data_origin" : {
"SET":"Domain_SET_ReferenceAlignment2024",
"KRAB":"Domain_KRAB_ReferenceAlignment2024",
"SSXRD":"Domain_SSXRD_ReferenceAlignment2024",
"ZF":"Domain_ZF_unit_ReferenceAlignment2024"
},
}

```


## Run the domain analysis on protein data

`uv run snakemake -s ../utils/process_stats_domain.smk  --jobs 1`

## Run the zf analysis on protein data

`uv run snakemake -s ../utils/process_zincfinger.smk  --jobs 1`

## Run the zf analysis on the dna sequence  coding for the protein data

`uv run snakemake -s ../utils/process_zincfinger_dna.smk  --jobs 1`

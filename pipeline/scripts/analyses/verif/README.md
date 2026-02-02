PROTEIN DOMAIN ANALYSIS ON PROTEOMES
====================================

## Settings 

We will use "uv" to run snakemake because it is the easisest way of managing python pacakges and versions especially when we hav no root permission.

Type :
```
uv init
uv add pandas
uv add hmmer
```

R has to be installed with the package Biostrings

## Configuration files


## Run the domain analysis on protein data

`uv run snakemake -s ../utils/process_stats_domain.smk  --jobs 1`

## Run the zf analysis on protein data

`uv run snakemake -s ../utils/process_zincfinger.smk  --jobs 1`

## Run the zf analysis on the dna sequence  coding for the protein data

`uv run snakemake -s ../utils/process_zincfinger_dna.smk  --jobs 1`

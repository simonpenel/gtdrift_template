Programs to be installed

- bbmap

- miniprot

Utilisation de busco

git clone de busco dans le home, puis

uv add ~/busco/


[project]
name = "analyse-busco-genomic"
version = "0.1.0"
description = "Add your description here"
readme = "README.md"
requires-python = ">=3.13"
dependencies = [
    "bio>=1.8.3",
    "busco",
    "hmmer>=3.4.0.2",
    "pandas>=3.0.5",
    "snakemake>=9.23.1",
]

[tool.uv.sources]
busco = { path = "../../../../../../../../beegfs/home/penel/busco" }


# Récupération des informations taxonomiques

## Génération du fichier ncbi_genome_assembly.txt
```
../utils/create_taxonomy_file.sh
```

Contenu du script:
```
#!/usr/bin/sh
echo -e "AssemblyName\tAssemblyAccession\tRefSeq_category\tSpeciesName\tSpeciesTaxid\tTaxid\tFtpPath_GenBank" > ncbi_genome_assembly.txt
esearch -db assembly -query eukaryota | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyName AssemblyAccession RefSeq_category SpeciesName SpeciesTaxid Taxid FtpPath_GenBank >> ncbi_genome_assembly.txt
sed 's/#//g' ncbi_genome_assembly.txt > tempfile && mv tempfile ncbi_genome_assembly.txt
```
## Génération du fichier ncbi_genome_assembly_taxon.txt
```
uv run ../utils/python/reformat_retrieve_taxonomy.py ncbi_genome_assembly.txt ncbi_genome_assembly_taxon.txt
```



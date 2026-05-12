# =======================================
# Join several  outputs. 

import pandas as pd
organisms_file = snakemake.input[0]
output_file = snakemake.output[0]
print(output_file)
print(organisms_file)
assembly=snakemake.params.accession
print(assembly)
df = pd.read_csv(organisms_file, sep='\t', header=0)
info = df.loc[df['Assembly Accession'] == assembly]
print(info)
info.to_csv(output_file, sep='\t',index=False)


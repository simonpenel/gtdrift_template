import json

# Function to load JSON files
def load_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

# Assign environment variables
globals().update(load_json("../environment_path.json"))


rule all:
    input:
        annoted="assemblies_annoted.json",
        not_annoted="assemblies_annoted.json"

rule split:
    """
    Split a list of assemblies  into annotated and not annotated 
    """
    priority: 1
    input:
        assemblies="assemblies.json",
        organisms=pathGTDriftData + "organisms_data"
    output:
        annoted="assemblies_annoted.json",
        not_annoted="assemblies_not_annoted.json"
    shell:
        """
        python3 ../utils/python/filter_assemblies.py -i {input.assemblies} -o {input.organisms}  -a {output.annoted} -n {output.not_annoted}
        """



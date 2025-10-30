import xml.etree.ElementTree as ET
import sys

# This script uses the ET library to read the xml file created as a result of the ncbi query.
# It creates a summary output text file named organisms_data and accessible in the data/resources
# directory.

tree = ET.parse(sys.argv[1])
root = tree.getroot()

meta_of_interest = ["non_chromosome_replicon_count","chromosome_count"]

with open(sys.argv[2], 'w') as writer:
    writer.write(f"Species Name\tTaxid\tAssembly Accession\tExisting Annotation\tExisting Protein Sequence\tURL")
    for meta in meta_of_interest:
        writer.write(f"\t"+meta)
    writer.write(f"\n")
    for document_summary in root.findall('.//DocumentSummary'):
        protein = False
        annotation = False
        source = 'GenBank'
        taxid = document_summary.find('Taxid').text
        meta = document_summary.find('Meta')
        dico = {} # dico categorie -> valeur
        stats = meta.find('Stats')  
        for stat in stats:
            for att in stat.attrib:
                categorie = stat.attrib[att]
                if categorie in meta_of_interest :
                    print(categorie+ " => "+stat.text)
                    if categorie in dico:
                        sys.exit("erreur")
                    dico[categorie] = stat.text
        accession = document_summary.find('AssemblyAccession').text
        if 'GCF' in accession:
            source = 'RefSeq'
        species_name = document_summary.find('SpeciesName').text

        properties = [elt.text for elt in document_summary.findall('PropertyList/string')]
        if any(prop in ['has_annotation', 'has_egap_annotation'] for prop in properties):
            protein = True
            annotation = True
        if document_summary.find(f"FtpPath_{source}") == None: # Pas de ftp pour les assemblages trop vieux ?
            url = None
        else:
            ftp = document_summary.find(f"FtpPath_{source}").text
            if ftp == None:
                url = None
            else:
                url = "https:" + ftp.strip().split(":")[1] + "/" + ftp.strip().split("/")[-1]
        status = document_summary.find('RefSeq_category').text 
        if status == 'representative genome' or  status == 'reference genome':
            line = f"{species_name}\t{taxid}\t{accession}\t{annotation}\t{protein}\t{url}"
            print(dico)
            for meta in meta_of_interest:
                val = dico[meta]
                line += f"\t{val}"
            line += f"\n"
            writer.write(line)

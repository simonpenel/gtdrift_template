import xml.etree.ElementTree as ET
import sys

# This script uses the ET library to read the xml file created as a result of the ncbi query.
# It creates a summary output text file named organisms_data and accessible in the data/resources
# directory.

tree = ET.parse(sys.argv[1])
flag_select_representative = True
if len(sys.argv) == 4 :
    if sys.argv[3] == "keep-all":
        flag_select_representative = False
        print("All accessions are selected")

root = tree.getroot()

with open(sys.argv[2], 'w') as writer:
    writer.write(f"Species Name\tTaxid\tAssembly Accession\tExisting Annotation\tExisting Protein Sequence\tURL\n")
    for document_summary in root.findall('.//DocumentSummary'):
        synonyms = document_summary.find('.//Synonym')
        protein = False
        annotation = False
        source = 'GenBank'
        taxid = document_summary.find('Taxid').text
        accession = document_summary.find('AssemblyAccession').text
        print("processing "+accession)
        if 'GCF' in accession:
            source = 'RefSeq'
        print("source = "+source)
        species_name = document_summary.find('SpeciesName').text

        properties = [elt.text for elt in document_summary.findall('PropertyList/string')]
        if any(prop in ['has_annotation', 'has_egap_annotation'] for prop in properties):
            protein = True
            annotation = True
        meta = document_summary.find('.//Meta')   
        if meta == None :
            ftp_sites = None
        else :
            ftp_sites = meta.find('.//FtpSites')
        if ftp_sites == None:
            print("No ftp for "+accession)
            url = None
        else :
            ftp_paths = ftp_sites.findall('.//FtpPath')
            ftp = None
            report = None
            for ftp_path in ftp_paths:
                print("ftp path = "+ftp_path.text)
                type = ftp_path.get("type")
                print("type = "+type)
                if type == source :
                    ftp = ftp_path.text
                if type == "Assembly_rpt" :
                    report = ftp_path.text
            if ftp == None:
                if report != None:
                    ftp = "/".join(report.split("/")[0:-1])

        if ftp == None:
            url = None
        else:
            url = "https:" + ftp.strip().split(":")[1] + "/" + ftp.strip().split("/")[-1]        
        # print("ftp path = "+ str(document_summary.find(f"FtpPath_{source}")))
        # if document_summary.find(f"FtpPath_{source}") == None: # Pas de ftp pour les assemblages trop vieux ?
        #     url = None
        # else:
        #     print("ftp path = "+ document_summary.find(f"FtpPath_{source}").text)
        #     ftp = document_summary.find(f"FtpPath_{source}").text
        #     if ftp == None:
        #         url = None
        #     else:
        #         url = "https:" + ftp.strip().split(":")[1] + "/" + ftp.strip().split("/")[-1]
        if flag_select_representative:
            status = document_summary.find('RefSeq_category').text 
            if status == 'representative genome' or  status == 'reference genome':
                line = f"{species_name}\t{taxid}\t{accession}\t{annotation}\t{protein}\t{url}\n"
                writer.write(line)
                print("ok "+accession+" "+status)
            else :
                print(accession+" "+status)
        else:
            line = f"{species_name}\t{taxid}\t{accession}\t{annotation}\t{protein}\t{url}\n"
            writer.write(line)
            # On ecrit aussi les synpnyme...
            # en effet  avec une requete 
            # esearch -db assembly -query   "GCA_947179515.1[Assembly Accession]" | efetch -format docsum
            # l'acession sera  GCF_947179515.1
            # il faut regadrer les synonyme pour pouvoir receuprer le GCA 
            # genbank = synonyms.find("Genbank")
            # if genbank != None:
            #     print("Genbank Accession "+genbank.text)

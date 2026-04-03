import argparse
import pandas
import subprocess

parser = argparse.ArgumentParser(description='Download proteins')

parser.add_argument('-d', '--data_organisms', type=str, required=True, help='organisms summary file')
parser.add_argument('-a', '--assembly', type=str, required=True, help='Assembly')
parser.add_argument('-o', '--output', type=str, required=True, help='File path')

args = parser.parse_args()

organisms_data = pandas.read_csv(args.data_organisms,sep="\t")
print(organisms_data)
selected = organisms_data[(organisms_data['Assembly Accession'] == args.assembly)]
print(selected)
url = list(selected["URL gff3 genome"])[0]
print(url)
subprocess.run(["wget",  url, "-O", args.output+".gz"])
subprocess.run(["gzip","-d", args.output+".gz"])
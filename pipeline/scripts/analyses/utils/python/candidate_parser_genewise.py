import pandas as pd
import argparse
from math import isnan

#pd.options.mode.copy_on_write = True
parser = argparse.ArgumentParser(description='Reads overview table in the csv format and returns best candidates for each locus')

parser.add_argument('-i', '--input', type=str, required=True, help='Overview table for prdm9')
parser.add_argument('-o', '--output', type=str, required=True, help='Processed file path')

args = parser.parse_args()

## Reading overview table for prdm9
table = pd.read_csv(args.input, sep=';', header=0)

## Sorting by chromosome then by start position
if len(table) > 0 :
    sorted = table.sort_values(by = ['Chromosome','Chr Start'])
else :
    sorted = table
## Getting rid of extra index column and resetting index
sorted.drop(columns='Unnamed: 0', inplace=True)
sorted.reset_index(drop=True, inplace=True)


## Columns to check if stop/frameshift in domain, and if truncated (stops in ZF)
sorted.insert(loc=len(sorted.columns)-3, column='Pseudogene (HMMER)', value='No')

## New columns for counts
sorted.insert(loc=len(sorted.columns)-2, column='Nb Stop/Frameshift', value=sorted['SET Stop/Frameshift']+sorted['KRAB Stop/Frameshift']+sorted['SSXRD Stop/Frameshift']+sorted['ZF Stop/Frameshift'])
sorted.insert(loc=len(sorted.columns)-2, column='Nb Introns', value=sorted['SET Intron']+sorted['KRAB Intron']+sorted['SSXRD Intron']+sorted['ZF Intron'])

## Columns to check if stop/frameshift in domain, and if truncated (stops in ZF)
#sorted.insert(loc=len(sorted.columns)-2, column='Pseudogene (HMMER)', value='No')
sorted.insert(loc=len(sorted.columns)-2, column='ZF Truncated', value='No')
for elt in sorted.index:
    if type(sorted['Stop/Shift Positions'][elt]) != float or not isnan(sorted['Stop/Shift Positions'][elt]):
        for pos in str(sorted['Stop/Shift Positions'][elt]).split(';'):
            if int(sorted['KRAB domain start'][elt]) < float(pos) < int(sorted['SET domain end'][elt]):
                sorted.loc[:,'Pseudogene (HMMER)'][elt] = 'Yes'
            if float(pos) < int(sorted['ZF domain start'][elt]) or int(sorted['ZF domain start'][elt]) < float(pos) < int(sorted['ZF domain end'][elt]):
                sorted.loc[:,'ZF Truncated'][elt] = 'Yes'

## Check if truncated (stop in ZF)

## Overlap function for comparison between two loci
def overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

## Fake Do-While loop: breaks once all comparisons have been made (k compared to length of index values)
k = 0
while True:

## Comparison between loci i and j
    i = 0
    while i < len(sorted.index)+2:
        j = 0
        
        while j < len(sorted.index)+2:
            if j == i or j not in sorted.index:
                j += 1
                continue
            if i not in sorted.index and i < len(sorted.index)+2:
                i += 1
                continue
            if i >= len(sorted.index)+2 or j >= len(sorted.index)+2:
                break
            
            itable = sorted.loc[i]
            jtable = sorted.loc[j]
            if overlap([int(itable['Chr Start']), int(itable['Chr End'])],
                    [int(jtable['Chr Start']), int(jtable['Chr End'])]) > 0 and itable['Chromosome'] == jtable['Chromosome'] and i != j:
                
                print("Overlap "+itable['SeqID'] + " with " + jtable['SeqID'])    
                
## First comparison: by score, each domain is attributed a score, if a locus has a better score the other one is dropped from table
                iscore = 0
                jscore = 0
                if itable['Nb SET domains'] > 0:
                    iscore += 10
                if itable['Nb KRAB domains'] > 0:
                    iscore += 10
                if itable['Nb SSXRD domains'] > 0:
                    iscore += 10
                if itable['Nb ZF domains'] > 0:
                    iscore += 1

                if jtable['Nb SET domains'] > 0:
                    jscore += 10
                if jtable['Nb KRAB domains'] > 0:
                    jscore += 10
                if jtable['Nb SSXRD domains'] > 0:
                    jscore += 10
                if jtable['Nb ZF domains'] > 0:
                    jscore += 1
                
                if iscore > jscore:
                    print(f"Dropping {jtable['SeqID']}: domains")
                    sorted.drop(index=j, inplace=True)
                    continue
                    
                elif jscore > iscore:
                    print(f"Dropping {itable['SeqID']}: domains")
                    sorted.drop(index=i, inplace=True)
                    i = j
                    j = 0
                    continue

## Second comparison: must account for truncation of domains
# ZF not accounted for: too variable

                itruncscore = 0
                jtruncscore = 0
                if itable['SET non-truncated'] > 0:
                    itruncscore += 10
                if itable['KRAB non-truncated'] > 0:
                    itruncscore += 10
                if itable['SSXRD non-truncated'] > 0:
                    itruncscore += 10

                if jtable['SET non-truncated'] > 0:
                    jtruncscore += 10
                if jtable['KRAB non-truncated'] > 0:
                    jtruncscore += 10
                if jtable['SSXRD non-truncated'] > 0:
                    jtruncscore += 10
                
                if itruncscore > jtruncscore:
                    print(f"Dropping {jtable['SeqID']}: truncated domains")
                    sorted.drop(index=j, inplace=True)
                    continue
                    
                elif jtruncscore > itruncscore:
                    print(f"Dropping {itable['SeqID']}: truncated domains")
                    sorted.drop(index=i, inplace=True)
                    i = j
                    j = 0
                    continue

## Third comparison: Number of stop codons, less stop codons is preferable (complete protein)
                elif itable['Nb Stop/Frameshift'] < jtable['Nb Stop/Frameshift']:
                    print(f"Dropping {jtable['SeqID']}: stops/frameshifts")
                    sorted.drop(index=j, inplace=True)
                    continue
                    
                elif jtable['Nb Stop/Frameshift'] < itable['Nb Stop/Frameshift']:
                    print(f"Dropping {itable['SeqID']}: stops/frameshifts")
                    sorted.drop(index=i, inplace=True)
                    i = j
                    j = 0
                    continue

## Fourth comparison: Pick a protein that is not a pseudogene
                elif itable['Pseudogene (HMMER)'] == 'No' and jtable['Pseudogene (HMMER)'] == 'Yes':
                    print(f"Dropping {jtable['SeqID']}: marked as pseudogene")
                    sorted.drop(index=j, inplace=True)
                    continue
                    
                elif jtable['Pseudogene (HMMER)'] == 'No' and itable['Pseudogene (HMMER)'] == 'Yes':
                    print(f"Dropping {itable['SeqID']}: marked as pseudogene")
                    sorted.drop(index=i, inplace=True)
                    i = j
                    j = 0
                    continue

## Fifth comparison: Protein length, longer proteins are more likely to be complete
                elif itable['Protein Length'] > jtable['Protein Length']:
                    print(f"Dropping {jtable['SeqID']} keeping {itable['SeqID']}: protein length")
                    sorted.drop(index=j, inplace=True)
                    continue
                    
                elif jtable['Protein Length'] > itable['Protein Length']:
                    print(f"Dropping {itable['SeqID']} keeping {jtable['SeqID']}: protein length")
                    sorted.drop(index=i, inplace=True)
                    i = j
                    j = 0
                    continue

## If passed, then both are identical: we remove one of them
                else:
                    print(f"Dropping {jtable['SeqID']} keeping {itable['SeqID']}: identity")
                    sorted.drop(index=j, inplace=True)
                    continue

            j += 1
        i += 1
    k += 1
    sorted.reset_index(drop=True, inplace=True)
    print(f'{k} / {len(sorted.index)}')
    if k >= len(sorted.index):
        break

# reorganise columns
# ------------------
sorted = sorted.astype({"Nb Stop/Frameshift": int, "Nb Introns": int, "SET Intron": int, "KRAB Intron": int, "SSXRD Intron": int, "ZF Intron": int })
sorted = sorted.astype({"SET Stop/Frameshift": int, "KRAB Stop/Frameshift": int, "SSXRD Stop/Frameshift": int, "ZF Stop/Frameshift": int })
sorted = sorted.astype({"SET non-truncated": int, "KRAB non-truncated": int, "SSXRD non-truncated": int, "ZF non-truncated": int })

sorted.to_csv(args.output, sep=';')

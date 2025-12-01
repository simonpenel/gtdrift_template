#!/bin/bash
# Download from NCBI protein sequences

# Paths
export genome_assembly=$1
export protOutput=$2
export symlink_directory=$3

echo "------------ ${genome_assembly} ------------"


if [[ ${genome_assembly} =~ GCA ]]; then
  PathLink=$(esearch -db assembly -query "${genome_assembly}[Assembly Accession]"    | efetch  -format docsum | xtract -pattern DocumentSummary -element FtpPath_GenBank)
elif [[ ${genome_assembly} =~ GCF ]]; then
  #echo "esearch -db assembly -query \"${genome_assembly}[Assembly Accession]\"" 
  PathLink=$(esearch -db assembly -query "${genome_assembly}[Assembly Accession]"   | efetch  -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq)
fi

PathLink=$(echo ${PathLink} | cut -d " " -f 1)

echo "------------ ${PathLink} ------------"

BASENAME=`basename ${PathLink}`

Path=https${PathLink:3}/${BASENAME}'_protein.faa.gz'

wget --spider ${Path}
if [[ $? == 0 ]]; then
wget ${Path} -O ${symlink_directory}${BASENAME}'_protein.faa.gz'
gzip -d ${symlink_directory}${BASENAME}'_protein.faa.gz'
cd ${symlink_directory}
ln -s ${BASENAME}'_protein.faa' ${protOutput}
else
echo "There is no protein data "
cd ${symlink_directory}
echo "THERE_IS_NO_PROTEIN_DATA" >  ${protOutput}
fi


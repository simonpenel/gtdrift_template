#!/bin/bash
# Download from NCBI cds sequences

# Paths
export genome_assembly=$1
export cdsOutput=$2
export symlink_directory=$3

echo "------------ ${genome_assembly} ------------"


if [[ ${genome_assembly} =~ GCA ]]; then
  PathLink=$(esearch -db assembly -query "${genome_assembly}[Assembly Accession]"   | efetch  -format docsum | xtract -pattern DocumentSummary -element FtpPath_GenBank)
elif [[ ${genome_assembly} =~ GCF ]]; then
  PathLink=$(esearch -db assembly -query "${genome_assembly}[Assembly Accession]"   | efetch  -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq)
fi

echo "------------ ${PathLink} ------------"

PathLink=$(echo ${PathLink} | cut -d " " -f 1)

BASENAME=`basename ${PathLink}`

Path=https${PathLink:3}/${BASENAME}'_cds_from_genomic.fna.gz'
wget --spider ${Path}

if [[ $? == 0 ]]; then

wget ${Path} -O ${symlink_directory}${BASENAME}'_cds_from_genomic.fna.gz'
gzip -d ${symlink_directory}${BASENAME}'_cds_from_genomic.fna.gz'

cd ${symlink_directory}
ln -s ${BASENAME}'_cds_from_genomic.fna' ${cdsOutput}
else
echo "There is no cds data "
cd ${symlink_directory}
echo "THERE_IS_NO_CDS_DATA" >  ${cdsOutput}
fi


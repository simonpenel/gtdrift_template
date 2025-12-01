#!/bin/bash
# Download from NCBI genome sequence

# Paths
export genome_assembly=$1
export genomeOutput=$2
export symlink_directory=$3
export storage=$4

echo "------------ ${genome_assembly} ------------"
echo "Storage : ${storage}"

if [[ ${genome_assembly} =~ GCA ]]; then
  PathLink=$(esearch -db assembly -query "${genome_assembly}[Assembly Accession]"   | efetch   -format docsum | xtract -pattern DocumentSummary -element FtpPath_GenBank)
elif [[ ${genome_assembly} =~ GCF ]]; then
  PathLink=$(esearch -db assembly -query "${genome_assembly}[Assembly Accession]"   | efetch    -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq)
fi

echo "------------initial pathlink  ${PathLink} ------------"
if [[ ${PathLink} == "" ]]; then
echo "Connection problem"
exit 0
fi

PathLink=$(echo ${PathLink} | cut -d " " -f 1)

echo "------------ ${PathLink} ------------"


BASENAME=`basename ${PathLink}`

Path=https${PathLink:3}/${BASENAME}'_genomic.fna.gz'

print
wget ${Path} -O ${symlink_directory}${BASENAME}'_genomic.fna.gz'
gzip -d ${symlink_directory}${BASENAME}'_genomic.fna.gz'

if [[ ${storage} == irods ]]; then
echo iput ${symlink_directory}${BASENAME}'_genomic.fna' /lbbeZone/home/penel/gtdrift/genome_seq/
iput ${symlink_directory}${BASENAME}'_genomic.fna' /lbbeZone/home/penel/gtdrift/genome_seq/
echo "Remove ${symlink_directory}${BASENAME}_genomic.fna"
rm ${symlink_directory}${BASENAME}'_genomic.fna'
fi


cd ${symlink_directory}
echo ${BASENAME}_genomic.fna > ${genomeOutput}
#if [[ ${storage} == irods ]]; then
#echo /lbbeZone/home/penel/gtdrift/genome_seq/${BASENAME}_genomic.fna > ${genomeOutput}
#else
#echo /lbbeZone/home/penel/gtdrift/genome_seq/${BASENAME}_genomic.fna > ${genomeOutput}
#ln -s ${BASENAME}'_genomic.fna' ${genomeOutput}
#fi

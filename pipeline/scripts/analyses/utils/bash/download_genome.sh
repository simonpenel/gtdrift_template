#!/bin/bash
# Download from NCBI genome sequence

# Paths
export organisms_data=$1
export genome_assembly=$2
export genomeOutput=$3
export symlink_directory=$4
export storage=$5

echo "------------ ${genome_assembly} ------------"
echo "Storage : ${storage}"


export Path=`grep $genome_assembly $organisms_data|cut -f6`
echo $Path
export nl=`echo $PATH|wc|awk '{print $1}'`
echo $nl
if [[ $nl -ne 1 ]]
then
echo "Duplication for $genome_assembly in  $organisms_data"
exit 1
fi

echo "[[[[ $Path]]]]"
export BASENAME=`basename ${Path}`
# #Path=https${PathLink:3}/${BASENAME}'_genomic.fna.gz'
export ftpPath=${Path}'_genomic.fna.gz'
echo $ftpPath
wget ${ftpPath} -O ${BASENAME}'_genomic.fna.gz'

gzip -d ${BASENAME}'_genomic.fna.gz'

if [[ ${storage} == irods ]]; then
echo iput -f {BASENAME}'_genomic.fna' /lbbeZone/home/penel/gtdrift/genome_seq/
iput -f  ${BASENAME}'_genomic.fna' /lbbeZone/home/penel/gtdrift/genome_seq/
echo "Remove ${BASENAME}_genomic.fna"
rm ${BASENAME}'_genomic.fna'
fi
# cd ${symlink_directory}
#echo ${BASENAME}_genomic.fna > ${symlink_directory}/${genomeOutput}

#cd ${symlink_directory}
if [[ ${storage} == irods ]]; then
echo "CREATE ${genomeOutput}"
echo /lbbeZone/home/penel/gtdrift/genome_seq/${BASENAME}_genomic.fna > ${genomeOutput}
else
mv ${BASENAME}'_genomic.fna' ${symlink_directory}
cd ${symlink_directory}
ln -s ${BASENAME}'_genomic.fna' ${genomeOutput}
fi

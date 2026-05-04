#!/usr/bin/sh
# select genomes with size < arg
echo $1
echo $2
head -1 $1 >  $1.selected.$2
head -1 $1 >  $1.excluded.$2
export acs=`cat $1 |cut -f3`
#echo $acs
for ac in $acs
do
#echo $ac
if [ -e ../../../../data_results_per_assembly/genome_assembly/$ac/genome_seq/genomic.fna.path ]
then
export path=`cat ../../../../data_results_per_assembly/genome_assembly/$ac/genome_seq/genomic.fna.path`
#echo $path
#ils -l  gtdrift/genome_seq/$path
#ils -l  gtdrift/genome_seq/$path|awk '{print $4}'
export size=`ils -l  gtdrift/genome_seq/$path|awk '{print $4}'`
echo "size = $size"
if [ $size -gt $2 ]
then
echo "superieur"
grep $ac $1 >> $1.excluded.$2
else
echo "inferieur"
grep $ac $1 >> $1.selected.$2
fi
else
echo "Erreur ../../../../data_results_per_assembly/genome_assembly/$ac/"
fi
done


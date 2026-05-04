#!/usr/bin/sh
# Test if softwares are available

# testing test

which programe_bidon
status=$?
echo "Test bidon : $status"
if  [ $status -eq 1 ]
then
echo "Ok"
else
echo "Error."
exit
fi



echo "MISSING SOFTWARE " > missing_soft 
echo "================ " >> missing_soft 



# Executables
# ===========

for soft in "genewise" "guix"  "formatdb" "blastp" "tblastn" "blastdbcmd" "python3" "pip3" "bidon"
do

which $soft
status=$?
echo "Test $soft : $status"
if  [ $status -eq 0 ]
then
echo "Ok"
else
echo "Error, $soft is missing"
echo $soft >> missing_soft 
fi
done


# Module phython
# ==============

which pip3
status=$?
if  [ $status -eq 0 ]
then
echo "pip3 Ok"
else
echo "Error, pip3 is missing"
exit
fi


echo "" >> missing_soft 
echo "MISSING MODULES " >> missing_soft 
echo "=============== " >> missing_soft 

for module in  "pandas" "bidon_mod" "ete3" "argparse"
do
echo test module $module
pip3 list|grep "^$module "
status=$?
echo "Module $module : $status"
if  [ $status -eq 0 ]
then
echo "Ok"
else
echo "Error, module $module is missing"
echo $module >> missing_soft 
fi
done






 



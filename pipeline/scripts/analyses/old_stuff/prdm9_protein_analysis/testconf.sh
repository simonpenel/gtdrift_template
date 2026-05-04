#!/usr/bin/sh
# Test if softwraes are available

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

for soft in "guix"  "formatdb" "blastp" "python3" "pip3" "bidon" "seqkit" "Rscript" "R"
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

for module in  "pandas" "bidon_mod" "ete3"
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

#module R
#biocmanager
#BiocManager::install("Biostrings")

for module in  "Biostrings"
do
echo test module $module
 Rscript -e  "library($module)"
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


 



#!/usr/bin/env bash

# uwemenzel@gmail.com  



## === Archive a GWAS project 




## +++ Calling:
#
# archive_gwas --id LIV_MULT3





## +++ Hardcoded settings & and defaults 

setfile=~/archive_settings.sh
if [ -s "${setfile}" ];then
  source ${setfile}  # command line paramters overwrite these settings  
else
  echo ""
  echo "  ERROR (archive_gwas.sh): Could not find the settings file \"${setfile}\"."
  echo ""
  exit 1  
fi


 
 
## +++ Command line parameters (override the settings in $setfile):

prog=$( basename "$0" )

if [ "$#" -lt 2 ]; then
  echo ""
  echo "  Usage: ${prog}"
  echo "         -i|--id <string>               no default"
  echo ""
  exit 1
fi
  

while [ "$#" -gt 0 ]
do
  case $1 in
      -i|--id)
          ident=$2
          shift
          ;;
      *)
          echo ""
	  echo "  Invalid argument: $1"
	  echo ""
	  exit 1
          ;;
  esac
  shift
done

 


## +++ Check if the variables are defined  

to_test=(ident)

for var in  ${to_test[*]}     
do
  if [ -z ${!var+x} ];then
    echo ""
    echo "  ERROR (archive_gwas.sh): mandatory variable $var is not defined."
    echo ""
  fi    
done




## +++ Check if folder exists:

if [ -d "${ident}" ]; then
  echo ""
  echo "  Going to compress the folder ${ident} ..." 
  echo ""
else 
  echo ""
  echo "  ERROR (archive_gwas.sh): Folder ${ident} not found."
  echo ""
  exit 1	  
fi


		

echo "  Keeping the following files outside:"
echo ""
nr_jma=$( ls ${ident}/${ident}_*.jma 2>/dev/null | wc -l )
if [ "${nr_jma}" -gt 0 ];then
  ls -1 ${ident}/${ident}_*.jma
  cp -i ${ident}/${ident}_*.jma .
else
  echo " No files containing pruned markers (.jma) found."
fi

nr_par=$( ls ${ident}/${ident}_gwas_params.txt 2>/dev/null | wc -l)
if [ "${nr_par=}" -gt 0 ];then
  ls -1 ${ident}/${ident}_gwas_params.txt
  cp -i ${ident}/${ident}_gwas_params.txt .
else
  echo " No parameter file found - this is probably an issue."
fi

nr_sig=$( ls ${ident}/${ident}_*gwas_signif.txt 2>/dev/null | wc -l )
if [ "${nr_sig}" -gt 0 ];then
  ls  -1 ${ident}/${ident}_*gwas_signif.txt
  cp -i ${ident}/${ident}_*gwas_signif.txt .
else
  echo " No file with significant markers found - this might be an issue."
fi
echo ""




## +++ Convert the time string for sbatch command below:

hours=$( expr $minutes / 60 )  
min=$( expr $minutes % 60 )    
if [ "$hours" -eq 0 ]; then
  time=${min}":00"
else
  if [ "${#min}" -eq 1 ]; then min="0"${min}; fi  
  time=${hours}":"${min}":00"    # this is the time for a single chromosome to run
fi  
echo "  Requested runtime for each chromosome: ${time}" | tee -a ${log}
echo "" | tee -a ${log}; echo "" | tee -a ${log}




## +++ Check if .tar.gz already exists

if [ -s "${ident}.tar.gz" ];then
  echo ""
  echo "  An archive ${ident}.tar.gz already exists."
  read -p "  Do you want to proceed anyway? (y/n):" -n 1 -r  
  echo
  if [[ $REPLY =~ ^[Yy]$ ]]; then 
    echo "  Starting compression"; echo
    rm -f ${ident}.tar.gz
  else
    echo; echo "  Bye."; echo
    exit 0
  fi
fi




## +++ Compress

account=$( echo $HOSTNAME | awk 'BEGIN{FS="-"} {print $1}' )
c_ident="ARCHIVE"
logf="${ident}_archive.log"

echo "  sbatch -A ${account} -p ${partition} -t ${time}  -J ${c_ident} -o ${logf} -e ${logf}  tar_gwas  --id ${ident}"

jobID=$(sbatch -A ${account} -p ${partition} -t ${time}  -J ${c_ident} -o ${logf} -e ${logf}  tar_gwas  --id ${ident} )  




## +++ Finish

echo ""
echo "  Done."
echo "" 






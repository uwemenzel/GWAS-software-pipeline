#!/usr/bin/env bash    



# uwemenzel@gmail.com



## === Remove samples from the genotype files but keep all variants 



## +++ Call:  
#
# remove_samples --genoid MF2 --new_genoid MF3 --sfile samples_to_remove.txt  





## +++ Hardcoded settings & and defaults 

setfile=~/remove_settings.sh
if [ -s "${setfile}" ];then
  source ${setfile}  # command line paramters overwrite these settings  
else
  echo ""
  echo "  ERROR (remove_samples.sh): Could not find the settings file \"${setfile}\"."
  echo ""
  exit 1  
fi


 
## +++ Command line parameters (override the settings in $setfile):

prog=$( basename "$0" )

if [ "$#" -lt 6 ]; then
  echo ""
  echo "  Usage: ${prog}"
  echo "         -g|--genoid <string>           no default"
  echo "         -n|--new_genoid                no default"
  echo "         -s|--sfile <file>              no default"     
  echo "         -c|--chr <int>[-<int>]         ${setfile}"  
  echo "         -m|--minutes <int>             ${setfile}"
  echo ""
  exit 1
fi

while [ "$#" -gt 0 ]
do
  case $1 in
      -g|--genoid)
          genoid=$2    
          shift
          ;;
      -n|--new_genoid)
          new_genoid=$2    
          shift
          ;;	  
      -s|--sfile)
          sfile=$2    
          shift
          ;;	    
      -c|--chr)
          chrom=$2
          shift
          ;;
      -m|--minutes)
          minutes=$2
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




## +++ Check if the variables are defined (including those defined in the settings file)    

to_test=(genoid new_genoid sfile chrom minutes partition minspace plink2_version genofolder)

for var in  ${to_test[*]}     
do
  if [ -z ${!var+x} ];then
    echo ""
    echo "  ERROR (${prog}): mandatory variable $var is not defined."
    echo ""
  fi    
done




## +++ Chromosomes  
      
cstart=$( echo $chrom | cut -d'-' -f 1 )
cstop=$( echo $chrom | cut -d'-' -f 2 )

if [[ ! ${cstart} =~ ^[0-9]+$ ]];then
  echo ""
  echo  "  ERROR (run_gwas.sh): Start chromosome is not valid: " ${cstart} 
  echo  "  			Correct syntax is e.g. --chrom 1-16"
  echo ""
  exit 1 
fi   

if [[ ! ${cstop} =~ ^[0-9]+$ ]];then
  echo ""
  echo  "  ERROR (run_gwas.sh): Stop chromosome is not valid: " ${cstop} 
  echo  "  			Correct syntax is e.g. --chrom 1-16"
  echo ""
  exit 1 
fi   

chromosomes=$( seq ${cstart} ${cstop} )



 
## +++ Header:

account=$( echo $HOSTNAME | awk 'BEGIN{FS="-"} {print $1}' ) 
log="${genoid}_${new_genoid}_remove.log"   
echo ""  > ${log}
echo ""   | tee -a ${log}
START=$(date +%s)      
echo -n "  "  | tee -a ${log}
date | tee -a ${log}
echo "  Account: ${account}" | tee -a ${log}
echo -n "  Operated by: " | tee -a ${log} 
whoami | tee -a ${log} 
echo "  Master logfile: ${log}" | tee -a ${log}
echo "" | tee -a ${log}
echo "  Running on chromosomes $cstart to $cstop" | tee -a ${log}  
echo "  Genotype input folder: ${genofolder}"  | tee -a ${log}
echo "  Genotype identifier: " ${genoid} | tee -a ${log}
echo "  New genotype identifier: " ${new_genoid} | tee -a ${log}
echo "  Samples to remove: ${sfile}"  | tee -a ${log}
echo "" | tee -a ${log}
echo "  Requested partition: ${partition}" | tee -a ${log}
echo "  Requested runtime per chromosome: ${minutes} minutes." | tee -a ${log}
echo "" | tee -a ${log}




## +++ Check available disk space:

echo -n "  Current working folder: " | tee -a ${log}
pwd | tee -a ${log}
space=$( df -k . | tail -1 | awk '{print $4}' )    
spac1=$( df -h . | tail -1 | awk '{print $4}' )  
echo "  Available disk space in this path: ${spac1}" | tee -a ${log}
echo "" | tee -a ${log}
if [ ${space} -lt ${minspace} ]; then   
    echo "" | tee -a ${log}
    echo "  Less than ${minspace} free disk space, consider using a different location." | tee -a ${log}
    echo "" | tee -a ${log}
    exit 1 
fi 




## +++ Check availability of input files and genotype files:

nr_found=0

pgentag=$( ls ${genofolder}/${genoid}*.pgen 2>/dev/null )
if [ "${pgentag}" != "" ];then
  nr_found=$((${nr_found} + 1))
  gtype="PGEN"
fi

bedtag=$( ls ${genofolder}/${genoid}*.bed 2>/dev/null )
if [ "${bedtag}" != "" ];then
  nr_found=$((${nr_found} + 1))
  gtype="BED"
fi

if [ "${nr_found}" -eq 0 ];then
  echo "" | tee -a ${log} 
  echo "  ERROR (${prog}): The folder ${genofolder} contains neither .pgen nor .bed genotype data." | tee -a ${log}
  echo "" | tee -a ${log} 
  exit 1  
fi

if [ "${nr_found}" -eq 2 ];then
  echo ""
  echo "  ERROR (${prog}): The folder ${genofolder} contains both .pgen and .bed genotype data." | tee -a ${log}  
  echo "                   Don't know what to do." | tee -a ${log} 
  echo "" | tee -a ${log} 
  exit 1
fi
      
echo "  Detected genotype: ${gtype}"  | tee -a ${log}
echo
 


if [ "${gtype}" == "PGEN" ];then

  for chrom in  ${chromosomes[*]}     
  do

    pgen_prefix="${genofolder}/${genoid}_chr$chrom"   

    psam=${pgen_prefix}".psam"  
    pvar=${pgen_prefix}".pvar"	
    pgen=${pgen_prefix}".pgen"	  

    if [ ! -f ${psam} ]; then
      echo "" | tee -a ${log}
      echo "  ERROR (${prog}): Input file '${psam}' not found." | tee -a ${log}
      echo "" | tee -a ${log}
      exit 1 
    fi  

    if [ ! -f ${pvar} ]; then
      echo "" | tee -a ${log}
      echo "  ERROR (${prog}): Input file '${pvar}' not found." | tee -a ${log}
      echo "" | tee -a ${log}
      exit 1 
    fi  

    if [ ! -f ${pgen} ]; then
      echo "" | tee -a ${log}
      echo "  ERROR (${prog}): Input file '${pgen}' not found." | tee -a ${log}
      echo "" | tee -a ${log}
      exit 1 
    fi	 

  done	

  echo "  All required genotype files (.pgen, .pvar, .psam) are available."  | tee -a ${log}
  
fi



if [ "${gtype}" == "BED" ];then

  for chrom in  ${chromosomes[*]}     
  do

    pgen_prefix="${genofolder}/${genoid}_chr$chrom"   

    bed=${pgen_prefix}".bed"  
    bim=${pgen_prefix}".bim"	
    fam=${pgen_prefix}".fam"	  

    if [ ! -f ${bed} ]; then
      echo "" | tee -a ${log}
      echo "  ERROR (${prog}): Input file '${bed}' not found." | tee -a ${log}
      echo "" | tee -a ${log}
      exit 1 
    fi  

    if [ ! -f ${bim} ]; then
      echo "" | tee -a ${log}
      echo "  ERROR (${prog}): Input file '${bim}' not found." | tee -a ${log}
      echo "" | tee -a ${log}
      exit 1 
    fi  

    if [ ! -f ${fam} ]; then
      echo "" | tee -a ${log}
      echo "  ERROR (${prog}): Input file '${fam}' not found." | tee -a ${log}
      echo "" | tee -a ${log}
      exit 1 
    fi	 

  done	

  echo "  All required genotype files (.bed, .bim, .fam) are available."  | tee -a ${log}
  
fi

echo "" | tee -a ${log}

if [ ! -f ${sfile} ]; then
  echo "" | tee -a ${log}   
  echo "  ERROR (${prog}): Input file '${sfile}' not found." | tee -a ${log}
  echo "" | tee -a ${log}
  exit 1 
else
  nr_entries=$( wc -l ${sfile} | awk '{print $1}' )
  echo "  The sample file has ${nr_entries} entries." | tee -a ${log}
fi  


# check sfile format (a little bit) 
nr_fields=$( awk '{print NF}' ${sfile} | sort | uniq -c | awk '{print $2}' )  
if [ "${nr_fields}" != "2" ];then
  echo "" | tee -a ${log}  
  echo "  ERROR (${prog}): Input file '${sfile}' must hav two columns: FID and IID." | tee -a ${log}
  echo "" | tee -a ${log}
  exit 1 
fi

echo "" | tee -a ${log}



## +++ Convert the time string for sbatch command below:

hours=$( expr $minutes / 60 )  
min=$( expr $minutes % 60 )    
if [ "$hours" -eq 0 ]; then
  time=${min}":00"
else
  if [ "${#min}" -eq 1 ]; then min="0"${min}; fi  
  time=${hours}":"${min}":00"    # time for a single chromosome to run
fi  







## +++ Run through chromosomes 

for chrom in  ${chromosomes[*]} 
do
     
  logchr="${genoid}_${new_genoid}_remove_chr${chrom}.log"   
  c_ident="RMOV-${chrom}"  		# unique jobID for each batch job
  			
  echo "sbatch -A ${account} -p ${partition} -t ${time} -J ${c_ident} -o ${logchr} -e ${logchr} \ "  | tee -a ${log}    
  echo "        remove_samples_chr --genoid ${genoid} --new_genoid ${new_genoid} --gtype ${gtype} --sfile ${sfile} --chr ${chrom}"  | tee -a ${log} 
  
  jobid=$( sbatch -A ${account} -p ${partition} -t ${time} -J ${c_ident} -o ${logchr} -e ${logchr} \
        remove_samples_chr --genoid ${genoid} --new_genoid ${new_genoid} --gtype ${gtype} --sfile ${sfile} --chr ${chrom} )

  echo "" | tee -a ${log}	  
   
done  


## +++ Finish  

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "" | tee -a ${log}
echo -n "  "  | tee -a ${log}
date | tee -a ${log}
echo "  Done." | tee -a ${log}
echo "" | tee -a ${log}
 









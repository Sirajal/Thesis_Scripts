#!/bin/bash
#$ -e /data4/siraj/logs
#$ -o /data4/siraj/logs
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -j y
#$ -pe openmpi 4
#$ -v PATH
#$ -v LD_LIBRARY_PATH
#$ -v PYTHONPATH

#Stringtie File
STRINGTIE="/data4/siraj/Stringtie/stringtie-1.3.3b.Linux_x86_64/stringtie"
#Get Filename
FILE_FULL=$(basename "$1")
FILE_NAME="${FILE_FULL%.*}"
#GTF Zebrafish
GTF="/data4/siraj/RNA_Seq_Data/Reference_Files/Danio_rerio.GRCz10.88_chr.gtf"
#Out Folder
OUT="/data4/siraj/RNA_Seq_Data/Stringtie_Out/"
#Make out folder
mkdir $OUT
#Run Stringtie
$STRINGTIE ${1} -p 4 -e -G ${GTF} -o ${DIR_OUT}/${FILE_NAME}.gtf

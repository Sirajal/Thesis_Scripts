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


#Kallisto File
KALLISTO="/data4/siraj/RNA_Seq_Data/Software/kallisto_linux-v0.43.1/kallisto"
#Kallisto Index
INDEX="/data4/siraj/RNA_Seq_Data/Indexes/Kallisto/"
#Out Folder
OUT="/data4/siraj/RNA_Seq_Data/Kallisto/"
#Make Out Folder
mkdir -p ${KALLISTO_OUT}/$1
#Run Kallisto
$KALLISTO quant -i ${INDEX}/GRCz10 -b 100 -t 4 -o ${OUT}/$1 ${2}/${1}"1.fastq.gz" ${2}/${1}"2.fastq.gz"

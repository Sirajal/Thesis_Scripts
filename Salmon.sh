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


#Salmon File
SALMON="/data4/siraj/RNA_Seq_Data/Software/Salmon-0.8.2_linux_x86_64/bin/salmon"
#Salmon Index
INDEX="/data4/siraj/RNA_Seq_Data/Indexes/Salmon/"
#Out Folder
OUT="/data4/siraj/RNA_Seq_Data/Salmon/"
#Make out folder
mkdir -p ${SALMON_OUT}/$1
#Run Salmon 
$SALMON quant -i ${INDEX}/GRCz10 -l A --seqBias --numBootstraps 100 -p 4 -o ${OUT}/$1 -1 ${2}/${1}"1.fastq.gz" -2 ${2}/${1}"2.fastq.gz"

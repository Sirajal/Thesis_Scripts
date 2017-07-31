#!/bin/bash
#$ -e /data4/siraj/logs
#$ -o /data4/siraj/logs
#$ -q all.q
#$ -cwd
#$ -S /bin/bash
#$ -j y
#$ -pe openmpi 1
#$ -v PATH
#$ -v LD_LIBRARY_PATH
#$ -v PYTHONPATH


#Modules to add
#module add bowtie/2-2.1.0
#module add tophat/2.0.10
#module add samtools

#Run Samtools Index
samtools index $1

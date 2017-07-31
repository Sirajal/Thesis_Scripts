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

#File name no extension
#https://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash
FILE_FULL=$(basename "$1")
FILE_NAME="${FILE_FULL%.*}"

#GTF File zebrafish
GTF="/data4/siraj/RNA_Seq_Data/Reference_Files/Danio_rerio.GRCz10.88_Chr.gtf"
#Output Folder
OUT="/data4/siraj/RNA_Seq_Data/HTSeq_Out/"
#Make output folder
mkdir $OUT

#Run BAM file to SAM file
#http://seqanswers.com/forums/showthread.php?t=24814
#Run HTSeq count
samtools view -h $1 | htseq-count -m intersection-nonempty -p yes -s no -f sam -r pos ${GTF} - ${OUT}/${FILE_NAME}.htseq.counts

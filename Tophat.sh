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

#Modules need
#module add bowtie/2-2.1.0
#module add tophat/2.0.10
#module add samtools

#GTF Zebrafish
GTF="/data4/siraj/RNA_Seq_Data/Reference_Files/Danio_rerio.GRCz10.88_chr.gtf"
#Out Folder
OUT="/data4/siraj/RNA_Seq_Data/Alignments"
#Bowtie Index
BT2_IDX="/data4/siraj/RNA_Seq_Data/Indexes/Bowtie2_Index/genome_ref"
#Transcript Index
TXI_IDX="/data4/siraj/RNA_Seq_Data/Indexes/Tophat_Transcript_Index/GRCz10/Danio_rerio.GRCz10.88_chr"
#Make out folder
mkdir -p ${OUT}/${1}
#Run Tophat
tophat2 -p 4 -o ${OUT}/${1} --b2-very-sensitive -G $GTF --transcriptome-index $TXI_IDX $BT2_IDX ${2}/${1}"1.fastq.gz" ${2}/${1}"2.fastq.gz"

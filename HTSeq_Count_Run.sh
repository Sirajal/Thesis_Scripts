#!/bin/bash

#BAM files are here
BAM_DIR="/data4/siraj/RNA_Seq_Data/Alignments/Complete_Alignments/"

#Loop BAM files
for i in $BAM_DIR/*.bam
do
	#HTSeq Count bam file
	qsub HTSeq_Counter.sh $i
done

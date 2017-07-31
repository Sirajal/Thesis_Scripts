#!/bin/bash

#BAM Files
BAM="/data4/siraj/RNA_Seq_Data/Alignments/Complete_Alignments/"

#Loop BAM files
for i in $BAM/*.bam
do
	#Run Samtools
	qsub Samtools.sh $i
done

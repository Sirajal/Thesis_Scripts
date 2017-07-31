#!/bin/bash

#BAM Files
BAM="/data4/siraj/RNA_Seq_Data/Alignments/Complete_Alignments/"

#Loop BAM Files
for i in $BAM/*.bam
do
	#Run Stringtie
	qsub Stringtie.sh $i
done

#!/bin/bash

#Bam Files
BAM="/data4/siraj/RNA_Seq_Data/raw_data/"
#Loop Bam Files
for i in $BAM/*_1.fastq.gz
do
	#Get Filename No Extension
	FILENAME=$(basename "$i")
	FILE_NOEXT="${FILENAME%.*}"
	FILE_NOEXT="${FILE_NOEXT%.*}"
	FILE_NOEXT=${FILE_NOEXT/%\_1/\_}
	#Run Tophat
	qsub Tophat.sh $FILE_NOEXT $BAM
done

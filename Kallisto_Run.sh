#!/bin/bash

#BAM Files
BAM="/data4/siraj/RNA_Seq_Data/raw_data/"

#Loop 1 files
for i in $BAM/*_1.fastq.gz
do
	#Get filename
	#https://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash
	FILENAME=$(basename "$i")
	FILE_NOEXT="${FILENAME%.*}"
	FILE_NOEXT="${FILE_NOEXT%.*}"
	FILE_NOEXT=${FILE_NOEXT/%\_1/\_}
	#Run Kallisto
	qsub Kallisto.sh $FILE_NOEXT $BAM
done

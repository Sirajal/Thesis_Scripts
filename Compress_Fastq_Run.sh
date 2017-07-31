#!/bin/bash

#Folder of fastq files
FQ_DIR="/data4/siraj/RNA_Seq_Data/raw_data/"

#Loop each fastq
for i in $FQ_DIR/*.fastq
do
	#Submit compress
	qsub Compress_Fastq.sh $i
done

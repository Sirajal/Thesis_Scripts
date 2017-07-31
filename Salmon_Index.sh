#!/bin/bash

#Zebrafish Transcripts
FASTA="/data4/siraj/RNA_Seq_Data/Reference_Files/GRCz10_Transcripts.fa"
#Salmon File
SALMON="/data4/siraj/RNA_Seq_Data/Software/Salmon-0.8.2_linux_x86_64/bin/salmon"
#Index File
INDEX="/data4/siraj/RNA_Seq_Data/Indexes/Salmon/"
#Make Index Folder
mkdir -p $INDEX
#Make Salmon Index
$SALMON index -i ${INDEX}/GRCz10 -t $FASTA

#!/bin/bash

#Transcript Fasta File
FASTA="/data4/siraj/RNA_Seq_Data/Reference_Files/GRCz10_Transcripts.fa"
#GTF Zebrafish
GTF="/data4/siraj/RNA_Seq_Data/Reference_Files/Danio_rerio.GRCz10.88_chr.gtf"
#Kallisto file
KALLISTO="/data4/siraj/RNA_Seq_Data/Software/kallisto_linux-v0.43.1/kallisto"
#Index File
INDEX="/data4/siraj/RNA_Seq_Data/Indexes/Kallisto/"
#Make output
mkdir -p $INDEX
#Run Kallisto
$KALLISTO index -i ${INDEX}/GRCz10 $FASTA

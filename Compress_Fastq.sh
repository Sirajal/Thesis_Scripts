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

#gz compress fastq file
gzip $1

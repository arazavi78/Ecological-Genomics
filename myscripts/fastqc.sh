#!/bin/bash


cd ~/Ecological-Genomics/myresults/fastqc

for file in /data/project_data/RS_ExomeSeq/fastq/edge_fastq/XPK*fastq.gz

do 

fastqc ${file} -o ./

done 





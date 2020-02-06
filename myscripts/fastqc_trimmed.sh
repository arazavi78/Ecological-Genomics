#!/bin/bash

mkdir fastqc_trimmed 

cd ~/Ecological-Genomics/myresults/fastqc_trimmed

for file in /data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/XPK*.cl.pd.fq

do 

fastqc ${file} -o ./

done 





#!/bin/bash   
 
cd /data/project_data/RS_ExomeSeq/fastq/edge_fastq  
 
for R1 in XPK*R1_fastq.gz  

do 
 
	R2=${R1/_R1_fastq.gz/_R2_fastq.gz}
	f=${R1/_R1_fastq.gz/}
	name=`basename ${f}`

	java -classpath /data/popgen/Trimmomatic-0.33/trimmomatic-0.33.jar org.usadellab.trimmomatic.TrimmomaticPE \
        -threads 1 \
# thread = Node 
        -phred33 \
# Note: the \ basically extends the line in bash 
         "$R1" \
         "$R2" \
         /data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/${name}_R1.cl.pd.fq \
         /data/project_data/RS_ExomeSeq/fastq/edge_fastq/unpairedcleanreads/${name}_R1.cl.un.fq \
         /data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/${name}_R2.cl.pd.fq \
         /data/project_data/RS_ExomeSeq/fastq/edge_fastq/unpairedcleanreads/${name}_R2.cl.un.fq \
        ILLUMINACLIP:/data/popgen/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:20 \
# this basically trimmes of anything with a phred score of lower than 20 (from the leading side and from the traiing side)
        TRAILING:20 \
        SLIDINGWINDOW:6:20 \ #didn't get that exactly 
        MINLEN:35 #minimum length 
 
done 


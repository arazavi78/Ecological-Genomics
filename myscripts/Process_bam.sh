#!/bin/bash

#This is where our output sam giles are going to get converted into binary format (bam)
#Then we are going to sort the bam files, remove the PCR duplicates, and index them.
#First, let's convert sam to bam and then we sory 

for f in ${output}/BWA/${mypop}*.sam


do 


  out=4{f/.sam/}
  sambamba-0.7.1-linux-static view -S --format=bam ${f} -o ${out}.bam
  samtools sort ${out}.bam -o ${out}.sorted.bam 

done 

#removing PCR replicates

for file in ${output}/BWA/${mypop}*.sorted.bam

do
#markdub deletes the PCR duplicates

  f=${file/.sorted.bam/}
  sambamba-0.7.1-linux-static markdub -r -t 1 ${file} ${f}}.sorted.rmdub.bam
  
  
done 



# to finish, we'll index our files.


for file in ${output}/BWA/${mypop}*.sorted.rmdub.bam


do 


  samtools index ${file}
  
  
done 









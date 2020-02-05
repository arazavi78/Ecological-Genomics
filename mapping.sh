#!/bin/bash


# path to the reference genome:

ref="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced"


#write a loop to map each individual within my population.

for forward in ${input}*R1.cl.pd.fq

do
  reverse=${forward/_R1.cl.pd.fq/_R2.cl.pd.fq}  
  # here we are replacing _R1.cl.pd.fq with _R2.cl.pd.fq
  f=${forward/_R1.cl.pd.fq/}
# here we do delete _R1.cl.pd.fq from the name of frowaard
  name= 'basename ${f}'
  #Basename only takes the end thiing in f
  # "'" means that you are giving a command within the variable.
  bwa mem -t 1 -M ${ref} ${forward} #{reverse} > ${output}/BWA/${name}.sam

done
#!/bin/bash


myrepo="/users/a/r/arazavi1/Ecological-Genomics/"

mypop="XPK"

input="/data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/${mypop}"

output="/data/project_data/RS_ExomeSeq/mapping/"

source ./mapping_final.sh

# Run the post-processing steps


source ./Process_bam.sh



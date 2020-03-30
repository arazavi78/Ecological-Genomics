#!/bash/bin


cd /data/project_data/RS_RNASeq/fastq/cleanreads

for file in NOR_02_H*R1.cl.fq
do 

salmon quant -i /data/project_data/RS_RNASeq/ReferenceTranscriptome/Pabies_HC_index -l A -r ${file} --validateMappings 
 --seqBias  -o /data/project_data/RS_RNASeq/salmon/HCmapping/${file}
done 





for file in XBM_07_C*R1.cl.fq

do 

salmon quant -i /data/project_data/RS_RNASeq/ReferenceTranscriptome/Pabies_HC_index -l A -r ${file} --validateMappings --seqBias  -o /data/project_data/RS_RNASeq/salmon/HCmapping/${file}

done


#!/bin/bash

# List of input files
files=(
    7C1_R1
    7C2_R1
    7C3_R1
    7C4_R1
    Biofilme1_R1
    Biofilme2_R1
    Biofilme3_R1
    Biofilme4_R1
    DMSO1_R1
    DMSO2_R1
    DMSO3_R1
    DMSO4_R1
    Planctonico1_R1
    Planctonico2_R1
    Planctonico3_R1
    Planctonico4_R1
    TC1_R1
    TC2_R1
    TC3_R1
    TC4_R1

)

# Loop through the files and run skewer
for file in "${files[@]}"; do
	docker run -v /mnt/d/Bioinfo/dados_completos_rna_seq:/home/ genomicpariscentre/samtools samtools view -bS -o /home/"$file"_mapped.bam /home/"$file"_mapped.sam
    docker run -v /mnt/d/Bioinfo/dados_completos_rna_seq:/home/ genomicpariscentre/samtools samtools sort /home/"$file"_mapped.bam -o /home/"$file"_mapped_sorted.bam
	docker run -v /mnt/d/Bioinfo/dados_completos_rna_seq:/home/ genomicpariscentre/samtools samtools index /home/"$file"_mapped_sorted.bam
done

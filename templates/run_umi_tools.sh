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
#docker run -it -v /mnt/d/Bioinfo/dados_completos_rna_seq:/home/ 83f49acd539b880cd258b27464d1480ffdee09f7df2935df25614e7f084d525e
# Loop through the files and run skewer
for file in "${files[@]}"; do
    umi_tools extract --stdin="$file".fastq.gz --extract-method=regex --bc-pattern=".{17,75}(?P<discard_1>AACTGTAGGCACCATCAAT)(?P<umi_1>.{12})(?P<discard_2>AGATCGGAAGAGCACACGTCT)(?P<discard_3>.*)" --log=processed.log --stdout "$file"_umi.fastq.gz
done

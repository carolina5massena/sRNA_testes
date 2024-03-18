#!/usr/bin/env nextflow
nextflow.enable.dsl=2


params.inputDir = "/mnt/d/Teste_Nextflow/arquivos_de_input"
params.outputDir = "/mnt/d/Teste_Nextflow"
params.genome = "/mnt/d/Teste_Nextflow/genoma_de_referencia/GCF_mixLupe.gtf" // valor padrão
params.index = "/mnt/d/Teste_Nextflow/genoma_de_referencia/index_GCF"


process umiToolsExtract {
    container "149243bc791fc6e729ac1daaedc04e2d8e4eb8f7e00b21dd0fd3482416ef53a3"

    input:
    path file

    output:
    path "*_umi.fastq"

    """
    umi_tools extract --stdin=${file} --extract-method=regex --bc-pattern='.{17,75}(?P<discard_1>AACTGTAGGCACCATCAAT)(?P<umi_1>.{12})(?P<discard_2>AGATCGGAAGAGCACACGTCT)(?P<discard_3>.*)' --log=processed.log --stdout ${params.outputDir}/umi_processados/${file}_umi.fastq.gz
    """
}

process bowtie {
    container "fd47fa53b2e"

    input:
    path file

    output:
    path "*_mapped.sam"

    """
    bowtie --threads 4 -v 2 -m 8 -a ${params.index} ${file} --sam ${params.outputDir}/bowtie_mapeados/${file}_mapped.sam
    """
}

process samtoolsView {
    container "genomicpariscentre/samtools"

    input:
    path file

    output:
    path "${file}_mapped.bam"

    """
    samtools view -bS -o ${params.outputDir}/samtools_convertidos/${file}_mapped.bam ${file}
    """
}

process samtoolsSort {
    container "genomicpariscentre/samtools"

    input:
    path file

    output:
    path "${file}_mapped_sorted.bam"

    """
    samtools sort ${file} -o ${params.outputDir}/samtools_ordenados/${file}_mapped_sorted.bam
    """
}

process samtoolsIndex {
    container "genomicpariscentre/samtools"

    input:
    path file

     output:
    path "${file}_mapped_sorted.bai"

    """
    samtools index ${file}
    """
}

process umiToolsDedup {
    container "149243bc791fc6e729ac1daaedc04e2d8e4eb8f7e00b21dd0fd3482416ef53a3"

    input:
    path file

    output:
    path "${file}_deduplicated.bam"

    """
    umi_tools dedup -I ${file} --output-stats=deduplicated -S ${params.outputDir}/umi_dedup_output/${file}_deduplicated.bam
    """
}

process featureCounts {
    container "5672d961627a"

    input:
    path file

    output:
    path "${file}_featurecounts.txt"

    """
    featureCounts -a ${params.genome} -O -g 'transcript_id' -o ${params.outputDir}/featureCounts/${file}_featurecounts.txt ${file}
    """
}

process fastqc_original {
    container "7b8f85bb68da"

    input:
    path file

    output:
    path "${file}.html"
    path "${file}.zip"

    """
    fastqc ${file} -o ${params.outputDir}/fastqc_output_original
    """
}
process fastqc_processado {  
    container "7b8f85bb68da"

    input:
    path file

    output:
    path "${file}.html"
    path "${file}.zip"

    """
    fastqc ${file} -o ${params.outputDir}/fastqc_output_processado
    """
}
workflow {

    sequencesFiles_ch = Channel.fromPath(params.inputDir + "/*").flatten()

    file("${params.outputDir}/featureCounts").mkdirs()
    file("${params.outputDir}/umi_dedup_output").mkdirs()
    file("${params.outputDir}/samtools_ordenados").mkdirs()
    file("${params.outputDir}/samtools_convertidos").mkdirs()
    file("${params.outputDir}/umi_processados").mkdirs()
    file("${params.outputDir}/bowtie_mapeados").mkdirs()
    file("${params.outputDir}/fastqc_output_processado").mkdirs()
    file("${params.outputDir}/fastqc_output_original").mkdirs()

    //fastqc_original(sequencesFiles_ch)
    umiToolsExtract(sequencesFiles_ch)
    //fastqc_processado(umiToolsExtract.out)
    bowtie(umiToolsExtract.out)
    samtoolsView(bowtie.out)
    samtoolsSort(samtoolsView.out)
    samtoolsIndex(samtoolsSort.out)
    umiToolsDedup(samtoolsSort.out)
    featureCounts(umiToolsDedup.out)
}
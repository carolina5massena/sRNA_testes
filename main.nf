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
    path "*_mapped.sam"

    """
    bowtie --threads 4 -v 2 -m 8 -a ${params.index} ${file} --sam ${params.outputDir}/bowtie_mapeados/${file}_mapped.sam
    """
}

process bowtie {
    container "fd47fa53b2e"

    input:
    path file

    output:
    path "${file}_mapped.sam"

    """
    bowtie --threads 4 -v 2 -m 8 -a ${params.index} ${file} --sam ${file}_mapped.sam
    """
}

process samtoolsView {
    container "genomicpariscentre/samtools"

    input:
    path file

    output:
    path "${file}_mapped.bam"

    """
    samtools view -bS -o ${file}_mapped.bam ${file}
    """
}

process samtoolsSort {
    container "genomicpariscentre/samtools"

    input:
    path file

    output:
    path "${file}_mapped_sorted.bam"

    """
    samtools sort ${file} -o ${file}_mapped_sorted.bam
    """
}

process samtoolsIndex {
    container "genomicpariscentre/samtools"

    input:
    path file

     output:
    path "${file}.bai"

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
    umi_tools dedup -I ${file} --output-stats=deduplicated -S ${file}_deduplicated.bam
    """
}

process featureCounts {
    container "5672d961627a"

    input:
    path file

    output:
    path "${file}_featurecounts.txt"

    """
    featureCounts -a ${params.genome} -O -g 'transcript_id' -o ${file}_featurecounts.txt ${file}
    """
}

process fastqc_original {
    container "7b8f85bb68da"

    input:
    path file

    """
    fastqc ${file} -o ${params.outputDir}/fastqc_output_original
    """
}
process fastqc_processado {
    container "7b8f85bb68da"

    input:
    path file

    """
    fastqc ${file} -o ${params.outputDir}/fastqc_output_processado
    """
}

workflow {
  
  
    file("${params.outputDir}/featureCounts").mkdirs()
    file("${params.outputDir}/umi_dedup_output").mkdirs()
    file("${params.outputDir}/samtools_ordenados").mkdirs()
    file("${params.outputDir}/samtools_convertidos").mkdirs()
    file("${params.outputDir}/umi_processados").mkdirs()
    file("${params.outputDir}/bowtie_mapeados").mkdirs()
    file("${params.outputDir}/fastqc_output_processado").mkdirs()
    file("${params.outputDir}/fastqc_output_original").mkdirs()



    // Create a channel with all files in params.inputDir
    sequencesFiles_ch = Channel.fromPath(params.inputDir + "/*").flatten()

    fastqc_original(sequencesFiles_ch)

    // Extract UMIs from the sequence files
    umiExtracted_ch = umiToolsExtract(sequencesFiles_ch)
    umiExtracted_ch.collectFile (
        storeDir: "${params.outputDir}/umi_processados"
    )

    // Align the sequences to the reference genome
    bowtieMapped_ch = bowtie(umiExtracted_ch)
    bowtieMapped_ch.collectFile (
        storeDir: "${params.outputDir}/bowtie_mapeados"
    )

    fastqc_processado(umiToolsExtract.out)

    // Convert SAM to BAM
    samView_ch = samtoolsView(bowtieMapped_ch)
    samView_ch.collectFile (
        storeDir: "${params.outputDir}/samtools_convertidos"
    )

    // Sort BAM files
    samSort_ch = samtoolsSort(samView_ch)
    samSort_ch.collectFile (
        storeDir: "${params.outputDir}/samtools_ordenados"
    )

    // Index sorted BAM files
    samIndex_ch = samtoolsIndex(samSort_ch)
    samIndex_ch.collectFile (
        storeDir: "${params.outputDir}/samtools_ordenados"
    )


    // Deduplicate BAM files
    umiDedup_ch = umiToolsDedup(samSort_ch)
    umiDedup_ch.collectFile (
        storeDir: "${params.outputDir}/umi_dedup_output"
    )

    // Count features
    featureCounts_ch = featureCounts(umiDedup_ch)
    featureCounts_ch.collectFile (
        storeDir: "${params.outputDir}/featureCounts"
    )

}

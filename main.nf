#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.inputDir = "inputs"
params.outputDir = "outputs"

process deduplicateUMIs {

    container "flomicsbiotech/umitools_last_version"

    input:
    path file

    output:
    path "deduplicated.bam"

    """
    umi_tools dedup -I ${file} --output-stats=deduplicated -S deduplicated.bam
    """
}

process getDataQualityReport {

    container "biocontainers/fastqc:v0.11.5"
    errorStrategy 'ignore'

    input:
    path file
    
    output:
    path "*.html"

    script:
    """
    fastqc -o . ${file}
    """
}


process extractUMIs {
    container "fredhutch/umi_tools:1.0.1"

    input:
    path file

    output:
    path "${file}_umi.fastq.gz"

    """
    umi_tools extract --stdin=${file} --extract-method=regex --bc-pattern='.{17,75}(?P<discard_1>AACTGTAGGCACCATCAAT)(?P<umi_1>.{12})(?P<discard_2>AGATCGGAAGAGCACACGTCT)(?P<discard_3>.*)' --log=processed.log --stdout ${file}_umi.fastq.gz
    """
}

workflow {
  
    // Create a channel with all files in params.inputDir
    sequencesFiles_ch = Channel.fromPath(params.inputDir + "/*").flatten()

    dataQualityReports_ch = getDataQualityReport(sequencesFiles_ch)
    dataQualityReports_ch.collectFile (
        storeDir: "${params.outputDir}/fastqc_reports"
    )

    extractedUMIs_ch = extractUMIs(sequencesFiles_ch)
    extractedUMIs_ch.collectFile (
        storeDir: "${params.outputDir}/extracted_umis"
    )
}
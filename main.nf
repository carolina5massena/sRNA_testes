params.inputDir = "inputs"
params.outputDir = "outputs"

process deduplicateUMIs {

    container = "flomicsbiotech/umitools_last_version"

    input:
    path file

    output:
    path "deduplicated.bam"

    """
    umi_tools dedup -I ${file} --output-stats=deduplicated -S deduplicated.bam
    """
}

workflow {
  
    // Create a channel with all files in params.inputDir
    umiSequencesFiles_ch = Channel.fromPath(params.inputDir + "/*")

    deduplicatedUMIs_ch = deduplicateUMIs(umiSequencesFiles_ch)

    deduplicatedUMIs_ch.collectFile (
        storeDir: "$params.outputDir"
    )

}
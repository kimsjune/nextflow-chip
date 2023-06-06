nextflow.enable.dsl=2

params.bed = "$projectDir/testCombine/*.bed"
params.bw = "$projectDir/testCombine/*.bw"

workflow {
    channel.fromFilePairs(params.bed, size:1)
    .combine(channel.fromFilePairs(params.bw, size:1))
    .view()
}

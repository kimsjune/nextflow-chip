nextflow.enable.dsl=2

params.index = "$projectDir/index/mm10.{1,2,3,4,rev.1,rev.2}.bt2"

workflow {
    channel.fromFilePairs(params.index, size:6)
    .view()
}

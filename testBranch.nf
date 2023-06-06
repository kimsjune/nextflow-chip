nextflow.enable.dsl=2

params.broadPeak = "$projectDir/testBroadPeak/*.broadPeak"

workflow {
    channel.fromFilePairs(params.broadPeak , size: 1)
        .branch { sample_id, broadPeak ->
            DMSO: sample_id.contains('DMSO')
                return [sample_id, broadPeak]
            nomatch: true
                return [sample_id, broadPeak]
        }
        .set{out}
    // out.DMSO.view()
    out.nomatch.view()
}
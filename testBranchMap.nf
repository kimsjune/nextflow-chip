nextflow.enable.dsl=2

params.input = "$projectDir/testBW/*.bw"
workflow {
    channel.fromFilePairs(params.input, size:1)
        .branch{sample_id, bw ->
            ChIP: sample_id.contains('ChIP')
                return [sample_id, bw]
        }
        .toSortedList( { a, b -> b[0] <=> a[0] } )
        .flatMap()
        .map {prefix, file ->
            def sample = prefix.toString().tokenize('_').get(0)
            tuple (sample, file) // The first string before the first underscore is the genotype. Need to find out how to make input come before ChIP for the heatmap.
        }
        .groupTuple()
        .view()
}
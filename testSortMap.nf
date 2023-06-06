nextflow.enable.dsl=2

params.sample = "$projectDir/testbam/*_{Rep1,Rep2}_sorted.bam"
input_ch = channel.fromFilePairs(params.sample, size : 1)

/*
workflow {
    input_ch
        .view()
}
*/
workflow {
    input_ch
    // https://nextflow-io.github.io/patterns/sort-filepairs-by-samplename/
    // explains how to sort and output the original tuple structure.
        .toSortedList( { a, b -> a[0] <=> b[0] } )
        .flatMap()
        .map {prefix, file ->
            def sample = prefix.toString().tokenize('_').get(0)
            def treatment = prefix.toString().tokenize('_').get(1)
            // def ip = prefix.toString().tokenize('_').get(2)
            def index = "${sample}_${treatment}"
            tuple(index, file)
            // thank god this finally works 
        }
        .groupTuple()
        //.view()
        .view()
}

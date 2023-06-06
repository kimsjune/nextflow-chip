nextflow.enable.dsl=2

params.fasta = "$projectDir/mm10.fa"
params.sample = "$projectDir/subsample/*_{R1,R2}.fastq"
params.index = "$projectDir/index/mm10.{1,2,3,4,rev.1,rev.2}.bt2"
params.flag = '0x2'
params.macsGenome = 'mm'
params.macsPeak = 'broad'
params.blacklist = "$projectDir/bed/mm10-blacklist.v2.bed"
params.genomeSize = 265283500
params.normMeth = 'RPGC'
params.repClass = "$projectDir/bed/mm10_rmsk_DNA_LowComplexity_Satellite_SimpleRepeat_sorted.bed" 
params.window = 1000
params.binSize = 38
params.heatmapCol = 'inferno'

//params.sampleChIP = "$projectDir/reads/*_{R1,R2}.fastq.gz"


sample_ch = channel.fromFilePairs(params.sample, size:2)
repClass_ch = channel.fromFilePairs(params.repClass, size:1)
index_ch = channel.fromFilePairs(params.index, size:6)
//sampleChIP_ch = channel.fromFilePairs(params.sampleChIP, size=2)


/* process BT2_INDEX {
    //memory 8.GB
    cpus 8

    input:
        tuple val(build), path(fasta)

    output:
        tuple val("$build"), path("$build*")
    script:
    """
    bowtie2-build --threads $task.cpus $fasta $build
    """

} */

process BT2_ALIGN {
    publishDir "$projectDir/sam", mode: "copy"
    cpus 8

    input:
        tuple val(sample), path(readPair)
        //tuple val(sampleChIP), path(chip_R1), path(chip_R2)
        tuple val(prefix), path(bt2) //path(bt2) is never actually used
    output:
        tuple val(sample), path("${sample}.sam"), emit:sam
        tuple val(sample), path("${sample}.log")
        //tuple val(sampleChIP), path("*.sam"), emit: chip
    script:
    """
    bowtie2 -t --sensitive-local -p $task.cpus -x $prefix -1 ${readPair[0]} -2 ${readPair[1]} -S ${sample}.sam 2> ${sample}.log
    """
}

process SAMTOOLS_VIEW {
    publishDir "$projectDir/bam", mode: "copy"
    cpus 4

    input:
        tuple val(sample), path(sam)
        //tuple val(sampleChIP), path(chipSam)
    output:
        tuple val(sample), path("*.bam"), emit:bam
        //tuple val(sampleChIP), path("*.bam"), emit: chip
    script:
    """
    samtools view -f $params.flag -b -@ $task.cpus $sam -o ${sample}.bam
    """
}

process SAMTOOLS_SORT {
    publishDir "$projectDir/bam", mode: "copy"
    cpus 4

    input:
        tuple val(sample), path(bam)
        //tuple val(sampleChIP), path(chipBam)
    output:
        tuple val(sample), path("${sample}_sorted.bam"), emit: bam
        //tuple val(sampleChIP), path("*_sorted.bam"), emit: chip
    script:
    """
    samtools sort -@ $task.cpus $bam -o ${sample}_sorted.bam
    """
}

process SAMTOOLS_INDEX {
    publishDir "$projectDir/bam", mode: "copy"
    cpus 4

    input:
        
        tuple val(sample), path(sortedBam)
        //tuple val(sampleChIP), path(chipSortedBam)
    output:
        tuple val(sample), path("${sample}*.bai"), emit: bai
        //tuple val(sampleChIP), path("*_sorted.bai")
    script:
    """
    samtools index $sortedBam 
    """
}
process BAMCOV {
    publishDir "$projectDir/bw", mode: "copy"
    cpus 8

    input:
        tuple val(sample), path(sortedBam), path(bai)

    output:
        tuple val(sample), path("*.bw")

    script:
    """
    bamCoverage -b $sortedBam -o ${sample}.bw -bs $params.binSize -bl $params.blacklist -p $task.cpus --effectiveGenomeSize $params.genomeSize --normalizeUsing $params.normMeth -e
    """
}

process MACS2 {
    publishDir "$projectDir/callpeak", mode: "copy"
    conda 'macs2.yml'

    input:
        tuple val(sample), path(lib)
                        // Or are these referred to by their indicies? 
   
    output:
        tuple val(sample), path("*_peaks.broadPeak"), emit: broadPeak
        tuple val(sample), path("*_peaks.gappedPeak")
        tuple val(sample), path("*_peaks.xls")
    script:
    """
    macs2 callpeak -t ${lib[0]} ${lib[2]} -c ${lib[1]} ${lib[3]} -f BAMPE -g $params.macsGenome -n $sample --$params.macsPeak
    """     
}

process BEDTOOLS_RM_BLK{
    //publishDir "$projectDir/callpeak", mode: "copy"

    input:
        tuple val(sample), path(broadPeak)
    output:
        tuple val(sample), path("*_rmblacklist.broadPeak")
    script:
    """
    bedtools intersect -a $broadPeak -b $params.blacklist -v >  ${sample}_rmblacklist.broadPeak
    """
}


process BEDTOOLS_INTERSECT {
    publishDir "$projectDir/bed", mode: "copy"
    input:
        tuple val(sample), path(broadPeak), val(repClass), path(bed) // only DMSO broadPeaks are used as inputs
        // need to combine these channels into one to get the Cartesian product of DMSO broadpeaks * repClass.bed files.
    output:
        tuple val("${sample}_${repClass}"), path("*.bed") //wow syntax important
        //https://github.com/nextflow-io/nextflow/issues/59
    script: // outputs DMSO peaks intersecting with repClass. Writes peaks.
    """
    bedtools intersect -wa -a $broadPeak -b $bed -u > ${sample}_${repClass}.bed
    """
}

process COMPMATRIX_ALL_PEAKS {
    input:
        tuple  val(region), path(broadPeak),val(sample), path(bw)
        
    output:
        tuple val(sample), path("*.npy.gz")
    script:
/*     bw input is sorted alphabetically. The order of bw input is the following: 
    ChIP DMSO rep1, input DMSO rep1,
    ChIP DMSO rep2, input DMSO rep2,
    ChIP GSK rep1, input GSK rep1
    ChIP GSK rep2, input GSK rep2
    There is no easy way of re-ordering this as I want from the workflow definition. Instead I am writing in the script here the order I want. */
    """
    computeMatrix scale-regions -R $broadPeak --sortRegions descend -p $task.cpus -S ${bw[1]} ${bw[3]} ${bw[5]} ${bw[7]} ${bw[0]} ${bw[2]} ${bw[4]} ${bw[6]} -o ${sample}.npy.gz --missingDataAsZero --skipZeros -b $params.window -a $params.window
    """
}

process HEATMAP_ALL_PEAKS{
    publishDir "$projectDir/svg", mode:"copy"
    input:
        tuple val(sample), path(matrix)

    output:
        tuple val(sample), path("*.svg")
    script:
    """
    plotHeatmap -z "" -m $matrix -o ${sample}.svg --colorMap $params.heatmapCol --heatmapHeight 14 --heatmapWidth 2 --startLabel "" --endLabel ""
    """
}

process  COMPMATRIX_PEAKSatRMSK {
    
    input:
        tuple val(genotype_repClass), path(repClassBed), val(genotype), path(bw)
         // only ChIP samples are used as inputs
                                      // Control Rep1, Rep2, Treatment Rep1, Rep2
    output:
        tuple val(genotype_repClass), path("*.npy.gz")

    script:
    """
    computeMatrix scale-regions -R $repClassBed --sortRegions descend -p $task.cpus -S ${bw[0]} ${bw[1]} ${bw[2]} ${bw[3]} -o ${genotype_repClass}.npy.gz --missingDataAsZero --skipZeros -b $params.window -a $params.window
    """

}

process PROFILE_PEAKSatRMSK {
    publishDir "$projectDir/svg", mode:"copy"
    input:
        tuple val(prefix), path(matrix)

    output:
        tuple val(prefix), path("*.svg")

    script:
    """
    plotProfile -z "" -m $matrix -o ${prefix}.svg --startLabel "" --endLabel "" --samplesLabel Rep1 Rep2 Rep1 Rep2 --colors black black blue blue --perGroup
    """

}



workflow {
    //BT2_INDEX(fasta_ch)
    BT2_ALIGN(sample_ch, index_ch.collect())

    SAMTOOLS_VIEW(BT2_ALIGN.out.sam)

    SAMTOOLS_SORT(SAMTOOLS_VIEW.out.bam)

    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

    SAMTOOLS_SORT.out.bam
        .combine(SAMTOOLS_INDEX.out.bai, by:0)
        .set{SAMTOOLS_SORTED_BAM_BAI}

    BAMCOV(SAMTOOLS_SORTED_BAM_BAI) // bw of all libraries 

    SAMTOOLS_SORT.out.bam // group bam files by ChIP and input libraries per genotype+treatment
        .toSortedList( { a, b -> a[0] <=> b[0] } )
        .flatMap()
        .map {prefix, file ->
            def genotype = prefix.toString().tokenize('_').get(0)
            def treatment = prefix.toString().tokenize('_').get(1)
            
            def fullPrefix = "${genotype}_${treatment}"
            tuple(fullPrefix, file)
            // thank god this finally works 
        }
        .groupTuple()
        .set{MACS_input_ch}

    MACS2(MACS_input_ch) // Biological replicates are pooled together

    BEDTOOLS_RM_BLK(MACS2.out.broadPeak) // genotype+treatment broadPeaks are filtered against blacklist bed file

    BEDTOOLS_RM_BLK.out // any treatment broadPeaks are not used to intersect with repClass
        .branch { sample_id, broadPeak ->
            DMSO: sample_id.contains('DMSO')
                return [sample_id, broadPeak]
            nomatch: true // this keeps the rest if needed.
                return [sample_id, broadPeak]
        }
        .set{BEDTOOLS_RM_BLK_BRANCH}

    BEDTOOLS_RM_BLK_BRANCH.DMSO
        .combine(repClass_ch)
        .set{broadPeak_repClass_comb}

    BEDTOOLS_INTERSECT(broadPeak_repClass_comb) 


    BAMCOV.out // bigwigs have to be grouped by input/ChIP and control/treatment libraries per genotype. Outputs a sorted list. Custom indices are declared in COMPMATRIX_PEAKS to reflect the order I want.
        .toSortedList( { a, b -> a[0] <=> b[0] } )
        .flatMap()
        .map {prefix, file ->
            def sample = prefix.toString().tokenize('_').get(0)
            tuple (sample, file) // The first string before the first underscore is the genotype. Need to find out how to make input come before ChIP for the heatmap.
        }
        .groupTuple()
        .set{COMPMATRIX_ALL_PEAKS_grouped}

    BEDTOOLS_RM_BLK_BRANCH.DMSO
        .combine(COMPMATRIX_ALL_PEAKS_grouped)
        .set{COMPMATRIX_ALL_PEAKS_input_ch}
        
    COMPMATRIX_ALL_PEAKS(COMPMATRIX_ALL_PEAKS_input_ch)

    HEATMAP_ALL_PEAKS(COMPMATRIX_ALL_PEAKS.out)

    BAMCOV.out 
        // first keep ChIP samples away from inputs
        .branch{sample_id, bw ->
            ChIP: sample_id.contains('ChIP')
                return [sample_id, bw]
        }
        // sort so that DMSO samples come before GSK (also works with Control and Treatment)
        .toSortedList( { a, b -> a[0] <=> b[0] } )
        .flatMap()

        // String before the first _ is the genotype that I want to group by
        .map {prefix, file ->
            def sample = prefix.toString().tokenize('_').get(0)
            tuple (sample, file)
        }
        .groupTuple()
        .set{BW_ChIP}

    BEDTOOLS_INTERSECT.out
        .combine(BW_ChIP)
        .set{COMPMATRIX_PEAKSatRMSK_input_ch}

    COMPMATRIX_PEAKSatRMSK(COMPMATRIX_PEAKSatRMSK_input_ch)

    PROFILE_PEAKSatRMSK(COMPMATRIX_PEAKSatRMSK.out) 
}
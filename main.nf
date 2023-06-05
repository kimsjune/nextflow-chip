nextflow.enable.dsl=2

params.fasta = "$projectDir/mm10.fa"
params.sample = "$projectDir/testreads/*_{R1,R2}.fastq.gz"
params.macsGenome = "mm"
params.macsPeak = "broad"
params.blacklist = "$projectDir/bed/mm10-blacklist.v2.bed"
params.genomeSize = "265283500"
params.normMeth = "RPGC"
params.repClass = "$projectDir/bed/mm10_rmsk_*_sorted.bed" 
params.window = "1000"
params.binSize = "38"

//params.sampleChIP = "$projectDir/reads/*_{R1,R2}.fastq.gz"


fasta_ch = channel.fromFilePairs(params.fasta, size:1)
sample_ch = channel.fromFilePairs(params.sample, size:2)
repClass_ch = channel.fromFilePairs(params.repClass, size:1)
//sampleChIP_ch = channel.fromFilePairs(params.sampleChIP, size=2)


process BT2INDEX {
    memory 8.GB

    input:
        tuple val(build), path(fasta)

    output:
        tuple val("$build"), path("$build*")
    script:
    """
    bowtie2-build $fasta $build
    """

}

process BT2ALIGN {
    // publishDir "$projectDir/sam"
    cpus 8

    input:
        tuple val(sample), path(readPair)
        //tuple val(sampleChIP), path(chip_R1), path(chip_R2)
        tuple val(prefix), path(bt2) //path(bt2) is never actually used
    output:
        tuple val(sample), path("*.sam")
        //tuple val(sampleChIP), path("*.sam"), emit: chip
    script:
    """
    bowtie2 -t --sensitive-local -p $task.cpus -x $prefix -1 ${readPair[0]} -2 ${readPair[1]} -S ${sample}.sam

    """
    /*bowtie2 -t --sensitive-local -p $task.cpus -x $prefix -1 $chip_R1 -2 $chip_R2 -S ${sampleChIP}.sam*/
}

process SAMTOOLS_VIEW {
    publishDir "$projectDir/aligned"
    cpus 4

    input:
        tuple val(sample), path(sam)
        //tuple val(sampleChIP), path(chipSam)
    output:
        tuple val(sample), path("*.bam")
        //tuple val(sampleChIP), path("*.bam"), emit: chip
    script:
    """
    samtools view -f 0x2 -b -@ $task.cpus $sam -o ${sample}.bam
    """
    //    samtools view -f 0x2 -b -@ $task.cpus $chipSam -o ${sampleChIP}.bam
}

process SAMTOOLS_SORT {
    publishDir "$projectDir/aligned", mode: "copy"
    cpus 4

    input:
        tuple val(sample), path(bam)
        //tuple val(sampleChIP), path(chipBam)
    output:
        tuple val(sample), path("*_sorted.bam")
        //tuple val(sampleChIP), path("*_sorted.bam"), emit: chip
    script:
    """
    samtools sort -@ $task.cpus $bam -o ${sample}_sorted.bam
    
    """
    //samtools sort -@ $task.cpus $chipBam -o ${sampleChIP}_sorted.bam
}

process SAMTOOLS_INDEX {
    publishDir "$projectDir/aligned", mode: "copy"
    cpus 4

    input:
        tuple val(sample), path(sortedBam)
        //tuple val(sampleChIP), path(chipSortedBam)
    output:
        tuple val(sample), path("*_sorted.bai")
        //tuple val(sampleChIP), path("*_sorted.bai")
    script:
    """
    samtools index $sortedBam

    """
    //samtools index $chipSortedBam
}
process BAMCOV {
    publishDir "$projectDir/bw", mode: "copy"
    cpus 8

    input:
        tuple val(sample), path(sortedBam)
    output:
        tuple val(sample), path("*.bw")

    script:
    """
    bamCoverage -b $sortedBam -o ${sample}.bw -bs $params.binSize -bl $params.blacklist -p $task.cpus --effectiveGenomeSize $params.genomeSize --normalizeUsing $params.normMeth -e
    """
}

process MACS2 {
    publishDir "$projectDir/callpeak", mode: "copy"

    input:
        tuple val(sample), path(lib)
                        // Or are these referred to by their indicies? 
   
    output:
        tuple val(sample), path("*_peaks.broadPeak"), emit: broadPeak
        tuple val(sample), path("*_peaks.gappedPeak")
        tuple val(sample), path("*_peaks.xls")
    script:
    """
    macs2 callpeak -t ${lib[0]} ${lib[1]} -c ${lib[2]} ${lib[3]} -f BAMPE -g $params.macsGenome -n $sample --$params.macsPeak
    """     
}

process BEDTOOLS_RM_BLK{
    publishDir "$projectDir/callpeak", mode: "copy"

    input:
        tuple val(sample), path(broadPeak)
    output:
        tuple val(sample), path("*_rmblacklist.broadPeak")
    script:
    """
    bedtools intersect -a $broadpeak -b $params.blacklist -v > \
    ${sample}_rmblacklist.broadPeak
    """
}


process BEDTOOLS_INTERSECT {
    publishDir "$projectDir/bed", mode: "copy"
    input:
        tuple val(sample), path(broadPeak) // only DMSO broadPeaks are used as inputs
        tuple val(repClass), path(bed)
    output:
        tuple val(sample_repClass), path("*.bed")
    script: // outputs DMSO peaks intersecting with repClass. Writes peaks.
    """
    bedtools intersect -b $bed -wa -a $broadPeak -u > ${sample}_${repClass}.bed
    """
}

process  COMPMATRIX_PEAKSatRMSK {
    
    input:
        tuple val(repClass), path(repClassBed)
        tuple val(genotype), path(bw) // only ChIP samples are used as inputs
                                      // Control Rep1, Rep2, Treatment Rep1, Rep2
    output:
        tuple val(genotype_repClass), path("*.npy.gz")

    script:
    """
    computeMatrix scale-regions -R $repClassBed --sortRegions descend -p $task.cpus -S ${bw[0]} ${bw[1]} ${bw[2]} ${bw[3]} -o $genotype_$repClass.npy.gz --missingDataAsZero --skipZeros -b $params.window -a $params.window
    """

}

process PROFILE_PEAKSatRMSK {

    input:
        tuple val(prefix), path(matrix)

    output:
        tuple val(prefix), path("*.svg")

    script:
    """
    plotProfile -z "" -m $matrix -o $prefix.svg --startLabel "" --endLabel "" --samplesLabel Rep1 Rep2 Rep1 Rep2 --colors black black blue blue --perGroup
    """

}




process COMPMATRIX_PEAKS {
    input:
        tuple val(sample), path(bw)
        tuple val(sample), path(broadPeak)
    output:
        tuple val(sample), path("*.npy.gz")
    script:
/*     bw input is sorted alphabetically. The order of bw input is the following: 
    ChIP GSK rep1, ChIP GSK rep2,
    ChIP DMSO rep1, ChIP DMSO rep2,
    input GSK rep2, input GSK rep2,
    input DMSO rep1, input DMSO rep2
    There is no easy way of re-ordering this as I want from the workflow definition. Instead I am writing in the script here the order I want. */
    """
    computeMatrix scale-regions -R $broadPeak -sortRegions descend -p $task.cpus -S ${bw[4]} ${bw[5]} ${bw[6]} ${bw[7]} ${bw[0]} ${bw[1]} ${bw[2]} ${bw[3]} -o ${sample}.npy.gz --missingDataAsZero --skipZeros -b $params.window -a $params.window
    """
}

workflow {
    BT2_INDEX(fasta_ch)
    BT2_ALIGN(sample_ch, BT2_INDEX.out.collect())
    SAMTOOLS_VIEW(BT2_ALIGN.out)
    SAMTOOLS_SORT(SAMTOOLS_VIEW.out)
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.)
    BAMCOV(SAMTOOLS_SORT.out) // bw of all libraries 

    // branch 1
    SAMTOOLS_SORT.out // group bam files by ChIP and input libraies per genotype+treatment
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
    BEDTOOLS_INTERSECT(BEDTOOLS_RM_BLK_BRANCH.DMSO, repClass_ch) 

    // branch 2
    BAMCOV.out // bigwigs have to be grouped by input/ChIP and control/treatment libraries per genotype. Outputs a sorted list. Custom indices are declared in COMPMATRIX_PEAKS to reflect the order I want.
        .toSortedList( { a, b -> a[0] <=> b[0] } )
        .flatMap()
        .map {prefix, file ->
            def sample = prefix.toString().tokenize('_').get(0)
            tuple (sample, file) // The first string before the first underscore is the genotype. Need to find out how to make input come before ChIP for the heatmap.
        }
        .groupTuple()
        .set{COMPMATRIX_PEAKS_input}

    COMPMATRIX_PEAKS(COMPMATRIX_PEAKS_input, BEDTOOLS_RM_BLK_BRANCH.DMSO)
    HEATMAP(COMPMATRIX_PEAKS.out)

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
    COMPMATRIX_PEAKSatRMSK(BEDTOOLS_INTERSECT.out, BW_ChIP )
    PROFILE_PEAKSatRMSK(COMPMATRIX_PEAKSatRMSK.out)


}
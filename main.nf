nextflow.enable.dsl=2

params.fasta = "$projectDir/dm6.fa"
params.sample = "$projectDir/testreads/*_{R1,R2}.fastq.gz"
params.macsGenome = "mm"
params.macsPeak = "broad"
params.blacklist = "$projectDir/bed/mm10-blacklist.v2.bed"
params.genomeSize = "265283500"
params.norm = "RPGC"
params.repClass = "$projectDir/bed/mm10_rmsk_*_sorted.bed"

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
        tuple val(sample), path(R1), path(R2)
        //tuple val(sampleChIP), path(chip_R1), path(chip_R2)
        tuple val(prefix), path(bt2) //path(bt2) is never actually used
    output:
        tuple val(sample), path("*.sam")
        //tuple val(sampleChIP), path("*.sam"), emit: chip
    script:
    """
    bowtie2 -t --sensitive-local -p $task.cpus -x $prefix -1 $R1 -2 $R2 -S ${sample}.sam

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

process MACS2 {
    publishDir "$projectDir/callpeak", mode: "copy"

    input:
        tuple val(sample), path(chipRep1), path(chipRep2)\
                           path(inputRep1), path(inputRep2) 
                        // Or are these referred to by their indicies? 
   
    output:
        tuple val(sample), path("*_peaks.broadPeak"), emit: broadPeak
        tuple val(sample), path("*_peaks.gappedPeak")
        tuple val(sample), path("*_peaks.xls")
    script:
    """
    macs2 callpeak -t $chipRep1 $chipRep2 -c $inputRep1 $inputRep2 -f BAMPE -g $params.macsGenome -n $sample --$params.macsPeak
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
        tuple val(sample), path(broadPeak)
        tuple val(repClass), path(bed)
    output:
        tuple val($sample_$repClass), path("*.bed")
    script:
    """
    bedtools intersect -b $bed -wa -a $broadPeak -u > ${sample}_${repClass}.bed
    """
}

process  COMPMATRIX_PEAKSatRMSK {

}

process PROFILE_PEAKSatRMSK {

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
    bamCoverage -b $sortedBam -o ${sample}.bw -bs 38 -bl $params.blacklist -p $task.cpus --effectiveGenomeSize $params.genomeSize --normalizeUsing RPGC -e
    """
}

process COMPMATRIX_PEAKS {
    input:
        tuple val(sample), path()
        tuple val(sample), path(broadPeak)
    output:
        tuple val(sample), path("*.npy.gz")
    script:
    """
    computeMatrix scale-regions -R $broadPeak -sortRegions descend -p $task.cpus -S $bw[0] $bw[0] $bw[0] $bw[0] -o ${sample}.npy.gz --missingDataAsZero --skipZeros -b 1000 -a 1000
    """
}

workflow {
    BT2_INDEX(fasta_ch)
    BT2_ALIGN(sample_ch, BT2_INDEX.out.collect())
    SAMTOOLS_VIEW(BT2_ALIGN.out)
    SAMTOOLS_SORT(SAMTOOLS_VIEW.out)
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.)

    // branch 1
    SAMTOOLS_SORT.out
        .toSortedList( { a, b -> a[0] <=> b[0] } )
        .flatMap()
        .map {prefix, file ->
            def sample = prefix.toString().tokenize('_').get(0)
            def treatment = prefix.toString().tokenize('_').get(1)
            def fullPrefix = "${sample}_${treatment}"
            tuple(fullPrefix, file)
            // thank god this finally works 
        }
        .groupTuple()
        .set{MACS_input_ch}
    MACS2(MACS_input_ch)
    BEDTOOLS_RM_BLK(MACS2.out.broadPeak)
    BEDTOOLS_INTERSECT(BEDTOOLS_RM_BLK.out, repClass_ch)
    COMPMATRIX_PEAKSatRMSK(BEDTOOLS_INTERSECT.out)
    PROFILE_PEAKSatRMSK(COMPMATRIS_PEAKSatRMSK.out)

    // branch 2
    BAMCOV(SAMTOOLS_SORT.out)
    BAMCOV.out
        .toSortedList( { a, b -> a[0] <=> b[0] } )
        // needs to be in reverse order... 
        .flatMap()
        .map {prefix, file ->
            def sample = prefix.toString().tokenize('_').get(0)
            tuple (sample, file)
        }
        .groupTuple()
        .set{COMPMATRIX_PEAKS_input}
    BEDTOOLS_RM_BLK.out
        .branch

        .set{BEDTOOLS_RM_BLK_DMSO}
    COMPMATRIX_PEAKS(COMPMATRIX_PEAKS_input, BEDTOOLS_RM_BLK_DMSO)
    HEATMAP(COMPMATRIX_PEAKS.out)

}
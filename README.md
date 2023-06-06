## Introduction
A Nextflow pipeline that demonstrates H3K27me3 ChIP-seq data analysis. 
Roughly speaking, the process is to:  
1. (Optional) Build a bowtie2 index. Or use pre-generated ones.
2. Align paired-end reads to the index.
3. Filter out concordantly mapped reads with samtools, then sort and index.
4. Identify broad peaks with MACS2.
5. Filter against blacklist.  


From here it diverges into two branches.  


First branch:  

1. Convert sorted bam files into bigWig using bamCoverage and normalized with RPGC<sup>*</sup>.
2. Compute matrix for vehicle treatment filtered peaks for all libraries per sample (biological replicates, input/ChIP, vehicle and drug treatment).
3. Plot a heatmap.  


Second branch:  


1. Intersect peaks with each RepeatMasker repClass.
2. Compute matrix for each repClass for all libraries of a sample except inputs.
3. Plot a profile.  

## To-do's and improvements
 - [ ] <sup>*</sup> for ChIP samples, `--scaleFactor` should be calcualted using `csaw` R package.
 - [x] Finish writing deeptools process and update workflow.
 - [ ] Include a new process that can generate RepeatMasker bed files, divided by repClass.
 - [x] Mermaid markdown.
 - [x] Change hard-coded variables into parameters whenever appropriate.
 - [x] Test out pipeline with a small genome or subsampled reads.
 - [ ] Separate main.nf into various modules.

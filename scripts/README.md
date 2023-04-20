# Folder Organization

## atac-seq
Contains scripts to:
* Align ATAC-seq reads, apply blacklists, deplete chrM reads
* Perform ATAC-seq shifts (if desired, as done by Buenrostro et al. (2013))
* Compute ATAC-seq peaks in 1 of 2 ways:
  1. Call ATAC-seq fragment pileup peaks (option `--broad` in MACS2)
  2. Call ATAC-seq cutsites (using `--shift` and `--extsize` in MACS2)
* Collapse ATAC-seq replicates with IDR
* Find enriched motifs in peaks using MEME-ChIP and TomTom
* R scripts to investigate expression changes proximal to differential changes in H3K27ac


## automateR
Contains source code for a package consisting of some simple custom functions that automates differential expression and differential enrichment analysis performed.

## chip-seq
Contains scripts to:
* Align ChIP-seq reads to human genome and apply necessary blacklists
* Check reproducibility of replicates by computing Pearson correlation coefficients with deepTools
* Call ChIP-seq peaks with MACS2
* Compute coverage density for a particular set of TSSs in a .bed file with deepTools
* R scripts for differential enrichment, GO term analysis

## hic
Contains scripts to:
* Remove reads mapping to chrY from HiC-Pro .matrix files (was necessary for my case in order for dcHiC to successfully complete SVD of interaction maps) and convert GRCh37 chromosome names to UCSC hg19 names
* Run dcHiC on interaction maps
* Perform gene expression/GO term analysis on genes annotated to changing windows

## integrated_analysis
Generate summary statistics for Wnt/pluripotency genes.

## preprocessing
Contains scripts to:
* Download a set of FASTQ files from a particular SRR accession
* Downsample FASTQ file(s) if desired
* Generate FastQC/MultiQC summary reports from FASTQ files
* Java script to parse BBDuk log files for the important statistics
* Trim reads with BBDuk

## rna-seq
Contains scripts to:
* Align RNA-seq data with HISAT2 and produce summary stats with RNASeQC
* Convert a set of .bam files to flat text count files with featureCounts
* Perform DGE analysis and produce summary MA, PCA, and volcano plots with DESeq2
* Other downstream analysis involving gene expression (e.g., heatmap visualization, GSEA)

## setup
Setup R environment, load all packages, get pluripotency/Wnt genes, get associated .bed file.
# Exploring Genome Reorganization and Gene Expression in Cancer Metastasis

### Author and Contact
Saurav Kiri, New College of Florida

For any inquiries, please contact [saurav.kiri19@ncf.edu](mailto:saurav.kiri19@ncf.edu).

## Overview
This GitHub repository contains all the code used my for my exploratory analysis thesis on comparing patterns of gene expression and epigenetic modifications in cancer metastasis as compared to primary tumors. The goal of my thesis was to identify possibly conserved pathways that are overexpressed or reprogrammed through structural re-organization after dissemination in an attempt to better understand the mechanims underlying metastasis. **Please see the folder `Scripts` for all code and descriptions of what is contained in each subfolder**.

In particular, I used the following GEO series datasets, including a mix of RNA-seq, ChIP-seq, ATAC-seq, and Hi-C sequencing:

[Integrative Hi-C maps of pancreatic cancer (GSE149103)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149103)

[Enhancer profiling in metastatic cancer (GSE98015)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98015)

[Specific chromatin landscapes and transcription factors couple breast cancer subtype with metastatic relapse to lung or brain (GSE129647)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129647)

I used a variety of pre-existing software tools to analyze ~ 3 TB of data in total downloaded from these repositories. In case you are interested in reproducing my pipeline, here are the steps I took for analyzing each dataset:

### RNA-seq
* Trim reads with BBDuk if necessary
* Align to GRCh37 with HISAT2
* Get gene counts with featureCounts
* Perform DGE analysis on each metastatic/primary tumor pairwise with DESeq2
* Run GSEA with clusterProfiler
  

### ChIP/ATAC-seq
* Trim reads with BBDuk if necessary
* Align to GRCh37 with BWA-MEM
* Remove low-quality reads, secondary alignments, unmapped reads, singletons (for PE data)
* Enforce hg19 blacklists
* Call peaks with MACS2
* Run differential enrichment analysis with DiffBind (`summits = 100` for ATAC-seq, `summits = 500` for H3K27ac ChIP-seq)
* Annotate regions to nearby genes
* Run GO term analysis on annotated genes

### Hi-C
* Process reads with HiC-Pro, removing low quality reads, singletons, and multi-mappers
* Create Hi-C interaction map at 40kb resolution
* Compute differential compartmentalization with dcHiC
* Annotate windows to genes, run GO term analysis on differentially compartmentalized genes

## Using this Code
This code is freely licensed under the the GNU Public License v3. A copy of this license is provided in this repository with NO WARRANTY. You are free to download, use, modify, and re-distribute all code in this GitHub. Just cite my page when you become famous.

## Installing automateR
To install my custom package `automateR`, simply download the `.tar.gz` file, open R to the directory containing the file, and type the following at the console:

``` r
install.packages("automateR_0.1.0.tar.gz", repos = NULL)
```

## Commandline dependencies
The following dependencies are necessary for properly running all semi-automated shell scripts:
* [SRA Toolkit](https://github.com/ncbi/sra-tools) (optional, if downloading data programmatically)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [Brian Bushnell's Decontamination Using K-mers (BBDuk)](https://sourceforge.net/projects/bbmap/)
* [HISAT2](https://github.com/DaehwanKimLab/hisat2)
* [featureCounts](https://subread.sourceforge.net)
* [BWA-MEM (part of BWA package)](https://github.com/lh3/bwa)
* [samtools](https://github.com/samtools/samtools)
* [bedtools](https://github.com/arq5x/bedtools2)
* [deepTools](https://github.com/deeptools/deepTools)
* [IDR](https://github.com/nboley/idr/tree/master)
* [MEME Suite](https://meme-suite.org/meme/doc/download.html)
* [MACS2](https://github.com/macs3-project/MACS)
* [HiC-Pro](https://github.com/nservant/HiC-Pro)
* [dcHiC](https://github.com/ay-lab/dcHiC)
project.dir <- "D:/SK/data"
ren.dir <- "ren-panc"
rna.dir <- "rna-seq/aligned"

setwd(file.path(project.dir, ren.dir, rna.dir))

source(file.path("C:", "Users", "jvons", "Documents", "NCF", "Thesis", "Scripts", "dge_analysis.r"))

featurecounts <- "ren_rna_rawcounts.txt"
cultures <- c("CAPAN1", "PANC1", "HPNE")

test.dge <- dge_analysis(featurecounts, 3, cultures)
test <- test.dge[[1]]

test.signif <- signifDE(test)

plotMA(test.dge[[1]])

volcano.plot(test.dge[[1]][, -ncol(test.dge[[1]])])

# The size factors range from 0.9592 - 1.08 - very good
    ## The fact that size factors are clustered indicates the sequencing of the libraries
    ## Was similiar in depth - our comparisons will be stronger
    ## Normalization technique:
    ### Step 1: for each gene, calculate the row-wise geometric means (i.e., geometric mean of its expression lvl across all samples)
    ### Step 2: For each gene in all samples, divide the raw count by its corresponding geometric mean
    ### Step 3: For each sample, take the median of normalized gene counts - these are the dfs
    ### Since the DESeq2 assumption is that most genes are not DE, the median normalized expression
    ### Should be around 1 by contruction

    # Pairwise comparison of all 3 cases
    ## Note that DESeq2's sample_1 vs. sample_2 indicates the LFC values are log2(sample_1/sample_2)
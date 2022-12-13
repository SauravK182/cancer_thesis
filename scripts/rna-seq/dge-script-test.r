featurecounts <- file.path(ren.data.dir, rna.data.dir, "ren_rna_rawcounts.txt")

coldata <- data.frame(Culture = rep(c("CAPAN1", "HPNE", "PANC1"), each = 3))
coldata$Culture <- factor(coldata$Culture, levels = c("CAPAN1", "PANC1", "HPNE"))

test.dge <- dge_analysis(featurecounts, coldata, contrasts = c("CAPAN1", "PANC1", "HPNE"))
test <- test.dge[[1]][[1]]
summary(test)

# See https://stackoverflow.com/questions/28543517/ for converting ENSEMBL to HGNC
require(EnsDb.Hsapiens.v79)
ensembl.genes <- rownames(test)
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79,
                             keys = ensembl.genes,
                             keytype = "GENEID",
                             columns = c("SYMBOL", "GENEID"))

# ENSEMBl IDs that map to NA HGNC symbols are excluded
# Need to manually match the HGNC symbol based on its assoc'd GENEID

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
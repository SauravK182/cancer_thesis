# Set up data frame for create.samplesheet()
tissue.cai <- rep(c("MDA-MB-231", "MDA-MB-231-BrM2", "MDA-MB-231-LM2"), each = 2)
factor.cai <- rep("h3k27ac", times = length(tissue.cai))
condition.cai <- rep(c("primary", "metastasis", "metastasis"), each = 2)
treatment.cai <- rep(NA, times = length(tissue.cai))
replicate.cai <- rep(c(1, 2), times = (length(tissue.cai) / 2))
sample.id.cai <- paste0(tissue.cai, replicate.cai)
cai.df <- data.frame(SampleID = sample.id.cai,
                     Tissue = tissue.cai,
                     Factor = factor.cai,
                     Condition = condition.cai,
                     Treatment = treatment.cai,
                     Replicate = replicate.cai)
# Set up directories
cai.bam.dir <- file.path(getwd(), cai.data.dir, "chip-seq/aligned")
cai.peaks.dir <- file.path(getwd(), cai.data.dir, "chip-seq/called_peaks")

# Call samplesheet to use for dba()
cai.samplesheet <- create_samplesheet(cai.df, cai.bam.dir, cai.peaks.dir)
cai.dba <- dba(sampleSheet = cai.samplesheet)

# Save peakset, explore heatmap
cai.peakset <- dba.peakset(cai.dba, bRetrieve = TRUE)
cai.cor <- dba.plotHeatmap(cai.dba)
cai.olap <- dba.overlap(cai.dba, mode = DBA_OLAP_RATE)

# Filter with blacklist and greylists
cai.filtered <- dba.blacklist(cai.dba, blacklist = FALSE, greylist = TRUE)
cai.filtered <- dba.blacklist(cai.filtered, blacklist = DBA_BLACKLIST_GRCH37, greylist = FALSE)

# Get summary of peak widths
summary(cai.filtered$binding[, 3] - cai.filtered$binding[, 2]) # min is 250, 1st Q is 823

# Count peaks
# cai.counts <- dba.count(cai.filtered, score = DBA_SCORE_RPKM, summits = 250)    # first start w/ 250
cai.counts <- dba.count(cai.filtered, score = DBA_SCORE_RPKM, summits = 500)    # now 500

cai.norm <- dba.normalize(cai.counts,
                          method = DBA_DESEQ2,
                          normalize = DBA_NORM_NATIVE)

# Exploratory plots of consensus set
dba.plotHeatmap(cai.norm)

# Set up contrast and analyze
cai.contrast <- dba.contrast(cai.norm,
                             design = ~ Tissue,
                             contrast = c("Tissue", "MDA-MB-231-LM2", "MDA-MB-231"))
cai.results <- dba.analyze(cai.contrast)
cai.DB <- dba.report(cai.results)
dba.plotHeatmap(cai.DB)
dba.plotMA(cai.results)

# Profile plot
cai.profile <- dba.plotProfile(cai.results)
dba.plotProfile(cai.profile)

#---------TESTING db_analysis() function------
test.analysis <- db_analysis(cai.samplesheet,
                             contrast.var = "Tissue",
                             summit.val = 500,
                             contrasts = list(c("MDA-MB-231-LM2", "MDA-MB-231")))
test.results.list <- test.analysis[[3]]
test.results <- test.results.list[[1]]
test.DB <- dba.report(test.results)

#-----------TEST RNA-SEQ---------------
# Set up coldata for DESeq2
cai.culture <- rep(c("BrM", "LM", "MB"), each = 3)
cai.type <- rep(c("MET", "MET", "PRT"), each = 3)
coldata.cai <- data.frame(Culture = cai.culture,
                          Type = cai.type)
coldata.cai$Culture <- factor(coldata.cai$Culture, levels = c("BrM", "LM", "MB"))
coldata.cai$Type <- factor(coldata.cai$Type, levels = c("MET", "PRT"))
coldata.cai$Condition <- factor(paste0(coldata.cai$Culture, coldata.cai$Type))

# Read in counts
cai.counts <- txt2counts(file.path(cai.data.dir, rna.data.dir, "cai_bc_rna_readcounts_raw.txt"))

# DDS
test.dds <- DESeqDataSetFromMatrix(countData = cai.counts,
                                   colData = coldata.cai,
                                   design = ~ Condition)
test.dseq <- DESeq(test.dds)
filter <- HTSFilter(test.dseq, s.len = 25, plot = FALSE)$filteredData
res <- results(filter, independentFiltering = TRUE, contrast = c("Condition", "BrMMET", "MBPRT"), alpha = 0.05)

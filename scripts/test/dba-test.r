# Set up columns of the sampleSheet object for DiffBind
tissue <- rep(c("HPNE", "PANC-1", "Capan-1"), each = 2)
sample.id <- paste0(tissue, rep(c(1, 2), times = 3))
factor <- rep("H3K27ac", times = length(sample.id))
condition <- rep(c("normal", "primary", "metastasis"), each = 2)
treatment <- rep(NA, times = length(sample.id))
replicate <- rep(c(1, 2), times = 3)

# Get bam reads
bam.dir <- file.path(main.data.dir, ren.data.dir, "chip-seq/aligned")
bam.files <- c()
for (type in unique(tissue)) {
    for (rep in c(1, 2)) {
    bam.file.name <- paste0(type, "_H3K27ac_ChIP-seq_", rep, ".bam")
    bam.file <- file.path(bam.dir, bam.file.name)
    bam.files <- c(bam.files, bam.file)
    }
}

# Get control files
control.id <- rep(c("HPNEc", "PANC-1c", "Capan-1c"), each = 2)
input.files <- c()
for (type in tissue) {
    input.name <- paste0(type, "_input.bam")
    input.loc <- file.path(bam.dir, input.name)
    input.files <- c(input.files, input.loc)
}

# Get peaks
peaks.dir <- file.path(main.data.dir, ren.data.dir, "chip-seq/called_peaks")
peaks <- c()
for (type in unique(tissue)) {
    for (rep in c(1, 2)) {
        peak.file.name <- paste0(type, "_H3K27ac_ChIP-seq_", rep, "_peaks.xls")
        peak.file.loc <- file.path(peaks.dir, peak.file.name)
        peaks <- c(peaks, peak.file.loc)
    }
}
peak.caller <- rep("macs", times = length(peaks))

# Assemble full sampleSheet
# Note, according to the following Bioconductor post, files must be in .bed format
# This means taking only first 6 columns and renaming to .bed: https://support.bioconductor.org/p/76216/
# Alternatively, the help file for dba (?dba) says to use the MACS .xls and have the PeakCaller column
# Be "macs" - this worked for me
ren.sample.sheet <- data.frame(SampleID = sample.id,
                               Tissue = tissue,
                               Factor = factor,
                               Condition = condition,
                               Treatment = treatment,
                               Replicate = replicate,
                               bamReads = bam.files,
                               ControlID = control.id,
                               bamControl = input.files,
                               Peaks = peaks,
                               PeakCaller = peak.caller)

ren.dba <- dba(sampleSheet = ren.sample.sheet)
dba.plotHeatmap(ren.dba)       # produces correlation heatmap from cross-correlation of each row of binding mat
# Can also use dba.plotPCA() to observe clustering of replicates
    # Use argument `label =` to add labels to the PCA

# MA plots of normalized and non-normalized data can be made as well with dba.plotMA

ren.peaks <- dba.peakset(ren.dba, bRetrieve = TRUE)

# Interestingly, for some reason, when specifying only blacklist, it uses the proper blacklist (GrCH37)
# But when I also specify greylist = TRUE, it states the genome build is BSgenome.Hsapiens.1000genomes.hs37d5 and fails to find a blacklist
# Must investigate further
# Building greylist here is probably redundant, since DiffBind does this when analyzing regardless
# The coverage removed from greylists is about identical to that done with DiffBind
blacklist.test <- dba.blacklist(ren.dba, blacklist = FALSE, greylist = TRUE)
blacklist.test2 <- dba.blacklist(blacklist.test, blacklist = DBA_BLACKLIST_GRCH37, greylist = FALSE)

# It seems that calling blacklist and greylists separately as above does the job quite well
# 75 consensus peaks were removed by enforcing both the blacklists and greylists

# One can use dba.overlap() to examine the number of consensus peaks overlapping in 1, 2,..., n samples
# Serves as a good QC plot as well - should be a shallow drop off, otherwise indicates lots of noise/outliers
# See the following two Bioconductor posts for choosing size for summits()
# https://support.bioconductor.org/p/100482/
# https://support.bioconductor.org/p/100170/
# https://www.biostars.org/p/9493721/ (for ATAC-seq)
# A recommendation is to examine the DBA object and look at the distribution of peak widths and choose e.g., the minimum or 1st quartile
# Of note, a recommendation for DiffBind is to use summits = 100 for ATAC-seq data, which
# will result in the creation of 201 bp peak windows. This is reasonable, since the found ATAC-seq fragment
# sizes by MACS is about 200 bp.

# Look at distribution of peak widths
summary(ren.dba$binding[, 3] - ren.dba$binding[, 2])
boxplot(ren.dba$binding[, 3] - ren.dba$binding[, 2])
# Here, the minimum is 199 bp and first quartile is 1583 bp - very diffuse peak widths

blacklist.test.count <- dba.count(blacklist.test2, summits = 500)   # Rory said he had previously used 500 for wide histone marks
blacklist.test.norm <- dba.normalize(blacklist.test.count,
                                     normalize = DBA_DESEQ2,
                                     library = DBA_NORM_NATIVE)
blacklist.norm <- dba.normalize(blacklist.test.norm, bRetrieve = TRUE)
# FrIPs are close to 10% - ENCODE recommends FRiP at least 1% (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/)

# For test blacklist + greylist
bltest.contrast <- dba.contrast(blacklist.test.norm,
                                 design = ~ Condition,
                                 contrast = c("Condition", "metastasis", "primary"))

# For test blacklist + greylist
bltest.analysis <- dba.analyze(bltest.contrast)
bltest.DB <- dba.report(bltest.analysis)        
# This time, 24550 peaks out of 39194 are DE - about 62.6% DE
# Note that for dba.report, one can pass a LFC threshold which will act similarly to DESeq2
# Namely, that the hypothesis test will be conducted for |LFC| > threshold, and p-values will be adjusted accordingly

# See https://support.bioconductor.org/p/79725/ for more

# Out of 39269 sites, 23847 were found to be differentially enriched for H3K27ac b/w PANC-1 and CAPAN-1
# Note however, this was just a preliminary analysis of going through the steps in the vignette, and
# I need to take more scrutiny when performing each individual step to better determine what it's doing
# Important note: CAPAN-1 and PANC-1 were isolated from different patients!
    # PANC-1 isolated from 56 yo white male (epithelioid carcinoma)
    # CAPAN-1 isolated from 40 yo white male (adenocarcinoma)

# MA plots post differential enrichment calling are also supported which will highlight intervals identified to be DE
# DiffBind also innately supports volcano plots with dba.plotVolcano()
# Additionally, the dba.plotHeatmap() also supports the results of DE enrichment via the `contrast` argument
    # to use only the intervals identified as DE
# See vignette also for concentration heatmaps

#----TEST FUNCTION----
ren.test <- db_analysis(ren.sample.sheet,
                        contrast.var = "Condition",
                        summit.val = 500,
                        contrasts = list(c("metastasis", "primary")))
ren.res <- ren.test[[3]][[1]]   # results seem to be the same as when I do it manually. n1

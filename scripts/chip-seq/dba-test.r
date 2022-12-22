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

# For now, use default when creating the dba.count object

# Next step: calculate binding matrix with scores based on read counts for every sample
ren.counts <- dba.count(ren.dba)
# FrIPs are close to 10% - ENCODE recommends FRiP at least 1% (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/)

# Now normalize with dba.normalize
ren.counts.norm <- dba.normalize(ren.counts, 
                                 method = DBA_DESEQ2, 
                                 normalize = DBA_NORM_NATIVE,
                                 library = DBA_LIBSIZE_FULL)
ren.norm <- dba.normalize(ren.counts.norm, bRetrieve = TRUE) # bRetrieve = TRUE will have normalization info returned

# Set up model design and contrast
ren.contrast <- dba.contrast(ren.counts.norm,
                             design = ~ Condition,
                             contrast = c("Condition", "primary", "metastasis"))

# Perform DE analysis - by default, runs DESeq2-based analysis
ren.analysis <- dba.analyze(ren.contrast)
dba.show(ren.analysis, bContrasts = TRUE)

# Out of 39269 sites, 23847 were found to be differentially enriched for H3K27ac b/w PANC-1 and CAPAN-1
# Note however, this was just a preliminary analysis of going through the steps in the vignette, and
# I need to take more scrutiny when performing each individual step to better determine what it's doing
# Important note: CAPAN-1 and PANC-1 were isolated from different patients!
    # PANC-1 isolated from 56 yo white male (epithelioid carcinoma)
    # CAPAN-1 isolated from 40 yo white male (adenocarcinoma)
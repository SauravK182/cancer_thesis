# Set up data frame for create.samplesheet()
tissue.rod <- rep(c("786-O", "786-M1A", "OS-RC2", "OS-LM1"), each = 2)
factor.rod <- rep("k27ac", times = length(tissue.rod))
condition.rod <- rep(rep(c("primary", "metastasis"), each = 2), times = 2)
treatment.rod <- rep(NA, times = length(tissue.rod))
replicate.rod <- rep(c(1, 2), times = (length(tissue.rod) / 2))
sample.id.rod <- paste0(tissue.rod, replicate.rod)
rod.df <- data.frame(SampleID = sample.id.rod,
                     Tissue = tissue.rod,
                     Factor = factor.rod,
                     Condition = condition.rod,
                     Treatment = treatment.rod,
                     Replicate = replicate.rod)
# Set up directories
rod.bam.dir <- file.path(getwd(), rod.data.dir, "chip-seq/aligned")
rod.peaks.dir <- file.path(getwd(), rod.data.dir, "chip-seq/called_peaks")

# Call samplesheet to use for dba()
rod.samplesheet <- create.samplesheet(rod.df, rod.bam.dir, rod.peaks.dir)
rod.dba <- dba(sampleSheet = rod.samplesheet)

# Save peakset
rod.peakset <- dba.peakset(rod.dba, bRetrieve = TRUE)
rod.cor <- dba.plotHeatmap(rod.dba)
rod.olap <- dba.overlap(rod.dba, mode = DBA_OLAP_RATE)

# Apply blacklist, greylist and count
rod.filtered <- dba.blacklist(rod.dba, blacklist = FALSE, greylist = TRUE)
rod.filtered <- dba.blacklist(rod.filtered, blacklist = DBA_BLACKLIST_GRCH37, greylist = FALSE)

# Get summary of peak widths, choose appropriate summit size
summary(rod.filtered$binding[, 3] - rod.filtered$binding[, 2])  # min = 199, 1st Q = 682
rod.counts <- dba.count(rod.filtered, summits = 200)

# Normalize, set up contrast and analyze
rod.norm <- dba.normalize(rod.counts,
                          method = DBA_DESEQ2,
                          normalize = DBA_NORM_NATIVE,
                          library = DBA_LIBSIZE_FULL)
rod.contrast <- dba.contrast(rod.norm,
                             design = ~ Tissue,
                             contrast = c("Tissue", "786-M1A", "786-O"))
rod.results <- dba.analyze(rod.contrast)
rod.DB <- dba.report(rod.results)
dba.plotMA(rod.results)

dba.plotBox(rod.results)    # boxplot of normalized read count in DE peaks

# Get correlation results from just the identified DE peaks
res.cor <- dba.plotHeatmap(rod.results)
# Strong correlation within cell lines, weak correlation across lines

########TEST ANNOTATE PEAKS########
require(ChIPpeakAnno)
require(EnsDb.Hsapiens.v75)
test.range <- toGRanges(EnsDb.Hsapiens.v75, feature = "gene")
test.anno <- annotatePeakInBatch(rod.DB,
                                 AnnotationData = test.range,
                                 output = "nearestBiDirectionalPromoter",
                                 bindingRegion = c(-5000, 5000))
# Need to figure out how to use GRCh37 coordinates/contig names and not hg19
# Get an error trying to annotatePeakInBatch since contigs are based in UCSC format ("chrN")
# Additionally, it seems the non-canonical contigs may or may not be an issue with annotation?
# Need to investigate further
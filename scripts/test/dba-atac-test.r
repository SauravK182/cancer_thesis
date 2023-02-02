# Set up df for ATAC-seq
atac.tis.ren <- rep(c("HPNE", "PANC-1", "Capan-1"), each = 2)
atac.factor.ren <- rep(c("ATAC-seq"), times = length(atac.tis.ren))
atac.cond.ren <- rep(c("normal", "primary", "metastasis"), each = 2)
atac.rep.ren <- rep(c(1, 2), times = length(atac.tis.ren) / 2)
atac.id.ren <- paste0(atac.tis.ren, atac.rep.ren)
atac.ren.df <- data.frame(SampleID = atac.id.ren,
                          Tissue = atac.tis.ren,
                          Factor = atac.factor.ren,
                          Condition = atac.cond.ren,
                          Replicate = atac.rep.ren)

# Set up reads and peaks dir
bam.dir <- file.path(main.data.dir, ren.data.dir, "atac-seq/aligned")
peaks.dir <- file.path(main.data.dir, ren.data.dir, "atac-seq/called_peaks")

# Get samplesheet
atac.ren.ss <- create_samplesheet(atac.ren.df, bam.dir, peaks.dir, atac = TRUE)

# Create DBA, perform blacklist
atac.ren.dba <- dba(sampleSheet = atac.ren.ss)
atac.ren.dba <- dba.blacklist(atac.ren.dba, blacklist = DBA_BLACKLIST_GRCH37, greylist = FALSE)

# Get peak widths
summary(atac.ren.dba$binding[, 3] - atac.ren.dba$binding[, 2])

# Get counts with summits = 100 (recommended for ATAC-seq)
#----------IMPORTANT NOTE---------
# DiffBind runs out of memory when trying to create peakset counts for the ATAC-seq data
# To account for this, going to try running the count function with DBA$config$yieldSize at a lower value
# DBA$config$yieldSize controls the number of reads processed at a single time, default 5000000
# Will try setting down to 1000000
atac.ren.dba$config$yieldSize <- 50000  # even at 50000, still getting memory errors
atac.ren.dba <-    dba.count(atac.ren.dba,  # will try using bUseSummarizeOverlaps = FALSE
                             summits = 100, # bUseSummarizeOverlaps = FALSE worked
                             score = DBA_SCORE_RPKM,
                             bUseSummarizeOverlaps = FALSE)

# Try now instead with summits = F to use native pileup
atac.ren.dba <-    dba.count(atac.ren.dba,  # will try using bUseSummarizeOverlaps = FALSE
                             summits = FALSE, # bUseSummarizeOverlaps = FALSE worked
                             score = DBA_SCORE_RPKM,
                             bUseSummarizeOverlaps = FALSE)

# Normalize
atac.ren.norm <- dba.normalize(atac.ren.dba,
                               method = DBA_DESEQ2,
                               normalize = DBA_NORM_RLE)
atac.ren.norm.info <- dba.normalize(atac.ren.norm, bRetrieve = TRUE)
atac.ren.norm.info

# Set up contrast and analyze
atac.ren.contrast <- dba.contrast(atac.ren.norm,
                                  design = ~ Condition,
                                  contrast = c("Condition", "metastasis", "primary"))
atac.ren.res <- dba.analyze(atac.ren.contrast)
atac.ren.DB <- dba.report(atac.ren.res)
atac.ren.DB     # 100806 ranges out of 143068 are differentially accessible with summits = 100
# With summits = FALSE, 102k ranges differentially accessible. will stick with summits = 100
dba.plotMA(atac.ren.res)


#----TEST FUN----
atac.ren.test <- db_analysis(atac.ren.ss,
                             contrast.var = "Condition",
                             summit.val = 100,
                             contrasts = list(c("metastasis", "primary")),
                             atac = TRUE)
atac.ren.res <- atac.ren.test[[3]][[1]] # 100806 ranges differentially accessible with summits = 100, same as above.

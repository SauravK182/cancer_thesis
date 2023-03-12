#-----------REN ATAC-SEQ---------
# Create minimal df for create_samplesheet()
tissue.ren <- rep(c("HPNE", "PANC-1", "Capan-1"), each = 2)
sample.id.ren <- paste0(tissue.ren, rep(c(1, 2), times = 3))
atac.factor.ren <- rep("ATAC-seq", times = length(sample.id.ren))
condition.ren <- rep(c("normal", "primary", "metastasis"), each = 2)
replicate.ren <- rep(c(1, 2), times = 3)
df.ren <- data.frame(SampleID = sample.id.ren,
                     Tissue = tissue.ren,
                     Factor = atac.factor.ren,
                     Condition = condition.ren,
                     Replicate = replicate.ren)

# Set up reads and peaks dir, create samplesheet
ren.atac.bam <- file.path(main.data.dir, ren.data.dir, "atac-seq/aligned")
ren.atac.peaks <- file.path(main.data.dir, ren.data.dir, "atac-seq/called_peaks")
atac.ren.ss <- create_samplesheet(df.ren, ren.atac.bam, ren.atac.peaks, atac = TRUE)

# Do DBA analysis with atac = TRUE
ren.atac <- db_analysis(atac.ren.ss,
                        contrast.var = "Condition",
                        summit.val = 100,
                        contrasts = list(c("metastasis", "primary")),
                        atac = TRUE)




#-----------CAI ATAC-SEQ---------
# Set up data frame for create_samplesheet()
tissue.cai <- rep(c("MDA-MB-231", "MDA-MB-231-BrM2", "MDA-MB-231-LM2"), each = 2)
atac.factor.cai <- rep("atac", times = length(tissue.cai))
condition.cai <- rep(c("primary", "metastasis", "metastasis"), each = 2)
replicate.cai <- rep(c(1, 2), times = (length(tissue.cai) / 2))
sample.id.cai <- paste0(tissue.cai, replicate.cai)
cai.df <- data.frame(SampleID = sample.id.cai,
                     Tissue = tissue.cai,
                     Factor = atac.factor.cai,
                     Condition = condition.cai,
                     Replicate = replicate.cai)

# Set up dirs, create samplesheet
cai.atac.bam <- file.path(main.data.dir, cai.data.dir, "atac-seq/aligned")
cai.atac.peaks <- file.path(main.data.dir, cai.data.dir, "atac-seq/called_peaks")
atac.cai.ss <- create_samplesheet(cai.df, cai.atac.bam, cai.atac.peaks, atac = TRUE)

# Do DBA analysis with atac = TRUE
contrasts.list <- list(c("MDA-MB-231-BrM2", "MDA-MB-231"),
                       c("MDA-MB-231-LM2", "MDA-MB-231"),
                       c("MDA-MB-231-BrM2", "MDA-MB-231-LM2"))
cai.atac <- db_analysis(atac.cai.ss,
                        contrast.var = "Tissue",
                        summit.val = 100,
                        contrasts = contrasts.list,
                        atac = TRUE)

#-------ANNOTATE PEAKS--------
ren.atac.anno <- anno_peak_gr37(dba.report(ren.atac[[3]][[1]]), output = "overlapping", region = NULL)
ren.atac.anno <- ren.atac.anno[!is.na(ren.atac.anno$feature), ]

cai.brm.atac.anno <- anno_peak_gr37(dba.report(cai.atac[[3]][[1]]), output = "overlapping", region = NULL)
cai.brm.atac.anno <- cai.brm.atac.anno[!is.na(cai.brm.atac.anno$feature), ]

cai.lm.atac.anno <- anno_peak_gr37(dba.report(cai.atac[[3]][[2]]), output = "overlapping", region = NULL)
cai.lm.atac.anno <- cai.lm.atac.anno[!is.na(cai.lm.atac.anno$feature), ]

anno.atac.list.full <- list(panc = ren.atac.anno,
                            brain_mb = cai.brm.atac.anno,
                            lung_mb = cai.lm.atac.anno)
#-----------REN H3K27AC---------
# Set up minimal df for samplesheet
tissue.ren <- rep(c("HPNE", "PANC-1", "Capan-1"), each = 2)
sample.id.ren <- paste0(tissue.ren, rep(c(1, 2), times = 3))
factor.ren <- rep("H3K27ac", times = length(sample.id.ren))
condition.ren <- rep(c("normal", "primary", "metastasis"), each = 2)
replicate.ren <- rep(c(1, 2), times = 3)
df.ren <- data.frame(SampleID = sample.id.ren,
                     Tissue = tissue.ren,
                     Factor = factor.ren,
                     Condition = condition.ren,
                     Replicate = replicate.ren)

# Create samplesheet
ren.bam <- file.path(main.data.dir, ren.data.dir, "chip-seq/aligned")
ren.peaks <- file.path(main.data.dir, ren.data.dir, "chip-seq/called_peaks")
ren.samplesheet <- create_samplesheet(df.ren, ren.bam, ren.peaks)

# Do DB analysis, isolate contrast
ren.chip <- db_analysis(ren.samplesheet,
                        contrast.var = "Condition",
                        summit.val = 500,
                        contrasts = list(c("metastasis", "primary")))
ren.chip.res <- ren.chip[[3]][[1]]




#---------RODRIGUES H3K27AC---------
# Set up data frame for create.samplesheet()
tissue.rod <- rep(c("786-O", "786-M1A", "OS-RC2", "OS-LM1"), each = 2)
factor.rod <- rep("k27ac", times = length(tissue.rod))
condition.rod <- rep(rep(c("primary", "metastasis"), each = 2), times = 2)
replicate.rod <- rep(c(1, 2), times = (length(tissue.rod) / 2))
sample.id.rod <- paste0(tissue.rod, replicate.rod)
rod.df <- data.frame(SampleID = sample.id.rod,
                     Tissue = tissue.rod,
                     Factor = factor.rod,
                     Condition = condition.rod,
                     Replicate = replicate.rod)

# Create samplesheet
rod.bam.dir <- file.path(main.data.dir, rod.data.dir, "chip-seq/aligned")
rod.peaks.dir <- file.path(main.data.dir, rod.data.dir, "chip-seq/called_peaks")
rod.samplesheet <- create_samplesheet(rod.df, rod.bam.dir, rod.peaks.dir)

# Do DB analysis
rod.chip <- db_analysis(rod.samplesheet,
                        contrast.var = "Tissue",
                        summit.val = 500,
                        contrasts = list(c("786-M1A", "786-O"), c("OS-LM1", "OS-RC2")))



#-------CAI H3K27AC--------------
# Set up data frame for create.samplesheet()
tissue.cai <- rep(c("MDA-MB-231", "MDA-MB-231-BrM2", "MDA-MB-231-LM2"), each = 2)
factor.cai <- rep("h3k27ac", times = length(tissue.cai))
condition.cai <- rep(c("primary", "metastasis", "metastasis"), each = 2)
replicate.cai <- rep(c(1, 2), times = (length(tissue.cai) / 2))
sample.id.cai <- paste0(tissue.cai, replicate.cai)
cai.df <- data.frame(SampleID = sample.id.cai,
                     Tissue = tissue.cai,
                     Factor = factor.cai,
                     Condition = condition.cai,
                     Replicate = replicate.cai)

# Create samplesheet
cai.bam.dir <- file.path(main.data.dir, cai.data.dir, "chip-seq/aligned")
cai.peaks.dir <- file.path(main.data.dir, cai.data.dir, "chip-seq/called_peaks")
cai.samplesheet <- create_samplesheet(cai.df, cai.bam.dir, cai.peaks.dir)

# Do DB analysis
contrasts.list <- list(c("MDA-MB-231-BrM2", "MDA-MB-231"),
                       c("MDA-MB-231-LM2", "MDA-MB-231"),
                       c("MDA-MB-231-BrM2", "MDA-MB-231-LM2"))
cai.chip <- db_analysis(cai.samplesheet,
                        contrast.var = "Tissue",
                        summit.val = 500,
                        contrasts = contrasts.list)


#-------ANNOTATE PEAKS--------
ren.chip.anno <- anno_peak_gr37(dba.report(ren.chip.res))
rod.786.chip.anno <- anno_peak_gr37(dba.report(rod.chip[[3]][[1]]))
rod.os.chip.anno <- anno_peak_gr37(dba.report(rod.chip[[3]][[2]]))
cai.brm.chip.anno <- anno_peak_gr37(dba.report(cai.chip[[3]][[1]]))
cai.lm.chip.anno <- anno_peak_gr37(dba.report(cai.chip[[3]][[2]]))

anno.chip.list.full <- list(panc = ren.chip.anno,
                       m1a_o = rod.786.chip.anno,
                       lm_rc = rod.os.chip.anno,
                       brain_mb = cai.brm.chip.anno,
                       lung_mb = cai.lm.chip.anno)

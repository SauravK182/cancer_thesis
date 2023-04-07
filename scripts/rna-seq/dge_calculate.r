#----------Ren et al., pancreatic cancer----------
# Get featurecounts file
featurecounts.ren <- file.path(ren.data.dir, rna.data.dir, "ren_rna_rawcounts.txt")

# Set up coldata for DESeq2
coldata.ren <- data.frame(Culture = rep(c("CAPAN1", "HPNE", "PANC1"), each = 3))
coldata.ren$Culture <- factor(coldata.ren$Culture, levels = c("CAPAN1", "PANC1", "HPNE"))

# Call DESeq2, isolate capan vs. panc comparison
ren.dge <- dge_analysis(featurecounts.ren, coldata.ren, contrasts = c("CAPAN1", "PANC1", "HPNE"))
capan.panc <- ren.dge[[1]][[1]]




#----------Rodrigues et al., ccRCC-----------------
# Get featurecounts txt file
featurecounts.rod <- file.path(rod.data.dir, rna.data.dir, "rod_rcc_rna_rawcounts.txt")

# Set up coldata for DESeq2
cell.lines <- rep(c("786",
                    "OS"), each = 4) %>% as.factor()
cancer.type <- rep(rep(c("MET", "PRT"), each = 2), 2)
coldata.rod <- data.frame(Culture = cell.lines, Type = cancer.type)
coldata.rod$Type <- factor(coldata.rod$Type, levels = c("MET", "PRT"))
coldata.rod$Condition <- factor(paste0(coldata.rod$Culture, coldata.rod$Type))

# Call DESeq2 for pairwise analysis of all conditions and isolate M1A vs. O and LM1 vs. RC2
# i.e., the metastasis vs. primary comparisons
dge.rod <- dge_analysis(featurecounts.rod, coldata.rod, formula.vec = c("Condition"))
rod.dge.list <- dge.rod[[1]]
m1a.o.comp <- rod.dge.list[[1]]
lm.rc.comp <- rod.dge.list[[6]]




#----------Cai et al., TN breast cancer-------------
# Get featurecounts file
featurecounts.cai <- file.path(cai.data.dir, rna.data.dir, "cai_bc_rna_readcounts_raw.txt")

# Set up coldata for DESeq2
cai.culture <- rep(c("BrM", "LM", "MB"), each = 3)
cai.type <- rep(c("MET", "MET", "PRT"), each = 3)
coldata.cai <- data.frame(Culture = cai.culture,
                          Type = cai.type)
coldata.cai$Culture <- factor(coldata.cai$Culture, levels = c("BrM", "LM", "MB"))
coldata.cai$Type <- factor(coldata.cai$Type, levels = c("MET", "PRT"))
coldata.cai$Condition <- factor(paste0(coldata.cai$Culture, coldata.cai$Type))

# Perform DESeq2 and isolate the three contrasts
dge.cai <- dge_analysis(featurecounts.cai, coldata.cai, contrast.var = "Culture", formula.vec = c("Culture"))
cai.dge.obj <- dge.cai[[2]]
dge.list.cai <- dge.cai[[1]]

brain.vs.primary <- dge.list.cai[[2]]
lung.vs.primary <- dge.list.cai[[3]]
brain.vs.lung <- dge.list.cai[[1]]




#--------Identifying genes DE across all samples-----
dge.list.full <- list(panc = capan.panc,
                      m1a_o = m1a.o.comp,
                      lm_rc = lm.rc.comp,
                      brain_mb = brain.vs.primary,
                      lung_mb = lung.vs.primary)

dge.list.signif <- lapply(dge.list.full, signifDE)
dge.list.signif.ens <- lapply(dge.list.signif, rownames)
dge.list.signif.hgnc <- lapply(dge.list.signif.ens, ensembl_to_gene)
signif.univ.ens <- purrr::reduce(dge.list.signif.ens, intersect)

# Get upregulated and downregulated genes for each experiment
dge.list.upreg <- lapply(dge.list.full, function(df) splitDE(df)[[1]])
dge.list.downreg <- lapply(dge.list.full, function(df) splitDE(df)[[2]])
dge.list.upreg.ens <- lapply(dge.list.upreg, rownames)
dge.list.downreg.ens <- lapply(dge.list.downreg, rownames)




#-------MAKING HEATMAP--------
full.df <- lapply(list(featurecounts.ren, featurecounts.rod, featurecounts.cai), txt2counts) %>%
            purrr::reduce(.f = function(df1, df2) {
                merge(df1, df2, by = 0, all = TRUE) %>%
                column_to_rownames(var = "Row.names")
            })
des.mat <- data.frame(Culture = c(rep(c("capan", "hpne", "panc"), each = 3),
                                  rep(c("m1a", "o", "lm1", "rc2"), each = 2),
                                  rep(c("mdamb", "lm", "br"), each = 3)),
                      Lab = c(rep("lab1", times = 9),
                              rep("lab2", times = 8),
                              rep("lab3", times = 9)))
dds <- DESeqDataSetFromMatrix(countData = full.df, colData = des.mat, design = ~ Culture)
full.dge <- DESeq(dds)
full.counts <- counts(full.dge, normalized = TRUE)

# Create Heatmap annotation
# See https://www.biostars.org/p/317349/
ann <- data.frame("Cell Line" = c(rep("Capan-1", 3),
                                  rep("HPNE", 3),
                                  rep("PANC-1", 3),
                                  rep("786-M1A", 2),
                                  rep("786-O", 2),
                                  rep("OS-LM1", 2),
                                  rep("OS-RC2", 2),
                                  rep("BrM2", 3),
                                  rep("LM2", 3),
                                  rep("MDA-MB", 3)),
                 "Morphology" = c(rep("Metastasis", 3),
                                  rep("Normal", 3),
                                  rep("Primary", 3),
                                  rep("Metastasis", 2),
                                  rep("Primary", 2),
                                  rep("Metastasis", 2),
                                  rep("Primary", 2),
                                  rep("Metastasis", 3),
                                  rep("Metastasis", 3),
                                  rep("Primary", 3)),
                 check.names = FALSE)

# Get colors for HeatmapAnnotation function
n <- length(unique(ann[, 1]))
# Code from https://stackoverflow.com/questions/15282580/
# require(RColorBrewer)
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# pie(rep(1,n), col=sample(col_vector, n))

# require(randomcoloR)
# palette <- distinctColorPalette(n)
# pie(rep(1, n), col = palette)

require(ComplexHeatmap)
col_vec <- viridis_pal(option = "plasma")(n)
names(col_vec) <- unique(ann[, 1])
morph_vec <- c("Metastasis" = "red", "Primary" = "#0026fd", "Normal" = "green")
col_vec <- list("Cell Line" = col_vec, "Morphology" = morph_vec)

colAnno <- HeatmapAnnotation(df = ann,
                             which = "column",
                             col = col_vec,
                             border = TRUE,
                             annotation_legend_param = list(at = c(unique(ann[, 1]), names(morph_vec))),
                             annotation_width = unit(c(1, 4), "cm"),
                             gap = unit(1, "mm"))

# z-score counts within sample and plot on Heatmap
heatmap.all <- t(apply(full.counts, 1, scale)) %>%
                na.omit() %>%
                Heatmap(cluster_rows = TRUE, cluster_columns = TRUE, show_column_names = FALSE,
                        name = "Z-score", show_row_names = FALSE,
                        top_annotation = colAnno)
save(heatmap.all, file = "C:/Users/jvons/Documents/NCF/Thesis/scripts/rna-seq/gene_all_heatmap.RData")

# Save heatmap
cairo_pdf("C:/Users/jvons/Documents/NCF/Thesis/Reports/gene_all_zacrossexp.pdf")
heatmap.all
dev.off()

# Take only genes DE in at least one comparison
de.genes <- lapply(dge.list.full, signifDE, lfc = 1) %>%
                lapply(function(deseq) rownames(as.data.frame(deseq))) %>%
                purrr::reduce(union)
heatmap.lfc1 <- t(apply(full.counts[de.genes, ], 1, scale)) %>%
                    na.omit() %>%
                    Heatmap(cluster_rows = TRUE, cluster_columns = TRUE, show_column_names = FALSE,
                            name = "Z-score", show_row_names = FALSE, show_row_dend = FALSE,
                            top_annotation = colAnno)

# Save heatmap
cairo_pdf("C:/Users/jvons/Documents/NCF/Thesis/Reports/gene_z_lfcthreshold_1.pdf")
heatmap.lfc1
dev.off()


# Create heatmap for LFC >= 2
de.genes <- lapply(dge.list.full, signifDE, lfc = 2) %>%
                lapply(function(deseq) rownames(as.data.frame(deseq))) %>%
                purrr::reduce(union)
heatmap.lfc2 <- t(apply(full.counts[de.genes, ], 1, scale)) %>%
                    na.omit() %>%
                    Heatmap(cluster_rows = TRUE, cluster_columns = TRUE, show_column_names = FALSE,
                            name = "Z-score", show_row_names = FALSE, show_row_dend = FALSE,
                            top_annotation = colAnno)

# Save heatmap
cairo_pdf("C:/Users/jvons/Documents/NCF/Thesis/Reports/gene_z_lfcthreshold_2.pdf")
heatmap.lfc2
dev.off()
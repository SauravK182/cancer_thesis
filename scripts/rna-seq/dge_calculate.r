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

# Get upregulated and downregulated genes for each experiment
dge.list.upreg <- lapply(dge.list.full, function(df) splitDE(df)[[1]])
dge.list.downreg <- lapply(dge.list.full, function(df) splitDE(df)[[2]])
dge.list.upreg.ens <- lapply(dge.list.upreg, rownames)
dge.list.downreg.ens <- lapply(dge.list.downreg, rownames)
features <- file.path(cai.data.dir, rna.data.dir, "cai_bc_rna_readcounts_raw.txt")
counts.df <- txt2counts(features)

cai.culture <- rep(c("BrM", "LM", "MB"), each = 3)
cai.type <- rep(c("MET", "MET", "PRT"), each = 3)
coldata.cai <- data.frame(Culture = cai.culture,
                          Type = cai.type)
coldata.cai$Culture <- factor(coldata.cai$Culture, levels = c("BrM", "LM", "MB"))
coldata.cai$Type <- factor(coldata.cai$Type, levels = c("MET", "PRT"))
coldata.cai$Condition <- factor(paste0(coldata.cai$Culture, coldata.cai$Type))

test.dge.cai <- dge_analysis(counts.df, coldata.cai, contrast.var = "Culture", formula.vec = c("Culture"), lfc = 0)
cai.dge.obj <- test.dge.cai[[2]]
dge.list.cai <- test.dge.cai[[1]]
brain.vs.primary <- dge.list.cai[[2]]
summary(brain.vs.primary)
lung.vs.primary <- dge.list.cai[[3]]

rld.cai <- rlog(cai.dge.obj)
plotPCA(rld.cai, intgroup = "Condition")

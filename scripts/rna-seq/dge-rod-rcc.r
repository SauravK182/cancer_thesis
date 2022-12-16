# Get featurecounts txt file
features <- file.path(rod.data.dir, rna.data.dir, "rod_rcc_rna_rawcounts.txt")

# Convert to counts df and set up colData for DESeq2
counts.df <- txt2counts(features)
cell.lines <- rep(c("786",
                    "OS"), each = 4) %>% as.factor()
cancer.type <- rep(rep(c("MET", "PRT"), each = 2), 2)
coldata <- data.frame(Culture = cell.lines, Type = cancer.type)
coldata$Type <- factor(coldata$Type, levels = c("MET", "PRT"))
coldata$Condition <- factor(paste0(coldata$Culture, coldata$Type))

# canc.dds <- DESeqDataSetFromMatrix(countData = counts.df,
#                                     colData = coldata,
#                                     design = ~ Condition)
# canc.dge <- DESeq(canc.dds)

# rld <- rlog(canc.dge, blind = TRUE)
# plotPCA(rld, intgroup = "Type")
# plotPCA(rld, intgroup = "Culture")

# There is relatively high within group variability for 

# dge.results <- lfcShrink(canc.dge, type = "normal", contrast = c("Type", "MET", "PRT"), alpha = 0.05, lfcThreshold = 0)
# plotMA(dge.results)

test.dge.rod <- dge_analysis(features, coldata, formula.vec = c("Condition"), lfc = 0)
# test.dge.rod <- dge_analysis(features, coldata, lfc = 0)
test.rod.list <- test.dge.rod[[1]]
m1a.o.comp <- test.rod.list[[1]]
summary(m1a.o.comp)
plotMA(m1a.o.comp)
m1a.o.signif <- rownames(signifDE(m1a.o.comp))

lm.rc.comp <- test.rod.list[[6]]
summary(lm.rc.comp)
plotMA(lm.rc.comp)
lm.rc.signif <- rownames(signifDE(lm.rc.comp))


#### Instead, subset data first and try ####
counts.df.subset <- counts.df[, 1:4]    # only gets data for 786 lines
coldata2 <- coldata[1:4, ]              # isolate rows for only 786
coldata2$Condition <- factor(coldata2$Condition, levels = c("786MET", "786PRT"))
test.rod.dge2 <- dge_analysis(counts.df.subset, coldata2, formula = c("Condition"), lfc = 0)

# Size factors are much better here than in the combined analysis
# 0.911 - 1.107 is the range of the 4 size factors
m1a.o.alone <- test.rod.dge2[[1]][[1]]
summary(m1a.o.alone)    # 223 upreg, 475 downreg at lfc = 1
                        # 4446 upreg, 4093 downreg at lfc = 0


# Find genes DE across all 3 cases
met.signif <- intersect(intersect(m1a.o.signif, lm.rc.signif), capan.panc.signif.ens)
length(met.signif)

# Create the DESeqDataSet
rod.dds <- DESeqDataSetFromMatrix(countData = counts.df,
                                  colData = coldata,
                                  design = test)

test.fun <- function(cell.lines, cancer.type = colnames(cell.lines)) {
    print(cancer.type)
}
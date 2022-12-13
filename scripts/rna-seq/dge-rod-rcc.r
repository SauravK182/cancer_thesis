# Get featurecounts txt file
features <- file.path(rod.data.dir, rna.data.dir, "rod_rcc_rna_rawcounts.txt")

# Convert to counts df and set up colData for DESeq2
counts.df <- txt2counts(features)
cell.lines <- rep(c("786-M1A",
                    "786-O",
                    "OS-LM1",
                    "OS-RC2"), each = 2) %>% as.factor()
cancer.type <- rep(rep(c("MET", "PRT"), each = 2), 2)
coldata <- data.frame(Culture = cell.lines, Type = cancer.type)
coldata$Type <- factor(coldata$Type, levels = c("PRT", "MET"))

test <- as.formula(paste0("~", "Culture"))

# Create the DESeqDataSet
rod.dds <- DESeqDataSetFromMatrix(countData = counts.df,
                                  colData = coldata,
                                  design = test)

test.fun <- function(cell.lines, cancer.type = colnames(cell.lines)) {
    print(cancer.type)
}
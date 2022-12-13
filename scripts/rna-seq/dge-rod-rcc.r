project.dir <- "D:/SK/data"
ren.dir <- "rod-rcc"
rna.dir <- "rna-seq/aligned"
setwd(file.path(project.dir, ren.dir, rna.dir))

source(file.path("C:", "Users", "jvons", "Documents", "NCF", "Thesis", "Scripts", "dge_analysis.r"))

counts.df <- txt2counts("rod_rcc_rna_rawcounts.txt")
cell.lines <- rep(c("786-M1A",
                    "786-O",
                    "OS-LM1",
                    "OS-RC2"), each = 2)
cancer.type <- rep(rep(c("MET", "PRT"), each = 2), 2)
coldata <- data.frame(Culture = cell.lines, Type = cancer.type)
coldata$Type <- factor(coldata$Type, levels = c("PRT", "MET"))
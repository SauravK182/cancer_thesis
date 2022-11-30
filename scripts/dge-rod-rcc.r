project.dir <- "D:/SK/data"
ren.dir <- "ren-panc"
rna.dir <- "rna-seq/aligned"
setwd(file.path(project.dir, ren.dir, rna.dir))

source(file.path("C:", "Users", "jvons", "Documents", "NCF", "Thesis", "Scripts", "dge_analysis.r"))
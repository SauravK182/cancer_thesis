# Read in diff comp files
get_comp <- function(name, getsubcomp = TRUE) {
    compdir <- file.path("E:/SK/data/", name, "dchic_res")
    comp <- file.path(compdir, "40kb/DifferentialResult/sample/fdr_result/")
    compfile <- file.path(comp, "differential.intra_sample_group.Filtered.pcQnm.bedGraph")
    subcompfile <- file.path(comp, "intra_sample_group.subcompartments.bedGraph")

    tryCatch(
        {
            message(paste0("Reading in from ", compfile, "..."))
            diffcomp <- read.delim(compfile)
            message("Done!")

            if (getsubcomp) {
                message(paste0("Reading in from ", subcompfile, "..."))
                subcomp <- read.delim(subcompfile)
                message("Done!")
                return(list(diffcomp = diffcomp, subcomp = subcomp))
            }
        },
        error = function(e) {
            stop("Differential compartment score file not found. Please try again")
        }
    )

    return(diffcomp)

}

# To assign either within or between compartment transitions
compswitch <- function(compdf, refcol = 5, expcol = 4) {
    aref <- compdf[, refcol] > 0
    aexp <- compdf[, expcol] > 0

    # Use mutate with case_when to add a column for compartment transition in vectorized format
    compdf.mut <- compdf %>%
                    mutate(transition = case_when(
                        aref & aexp ~ "Within A Transition",
                        aref & !aexp ~ "A to B Transition",
                        !aref & aexp ~ "B to A Transition",
                        !aref & !aexp ~ "Within B Transition"
                    ))
    
    return(compdf.mut)
}

panc.comp <- get_comp("ren-panc")
rcc.comp <- get_comp("rod-rcc", getsubcomp = FALSE)

# Assign transition states to compartment scores
# Note diff compartmentalization is called at FDR 0.10 threshold
hic.list <- list(panc = panc.comp[["diffcomp"]], m1a_o = rcc.comp)
hic.gr.list <- lapply(hic.list, function(df) {
    df %>%
        compswitch() %>%
        (function(x) {
            for (i in seq_len(nrow(x))) {
                x$chr[i] <- gsub(pattern = "chr", replacement = "", x = x$chr[i]) # convert to GRCh37 chr names
            }
            return(x)
        }) %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE)
})

# Annotate diff comp regions, reporting all genes overlapping each region
hic.anno.list <- lapply(hic.gr.list, function(gr) {
    gr.anno <- anno_peak_gr37(gr, output = "overlapping", region = NULL, select = "all")
    gr.omit <- gr.anno[!is.na(gr.anno$feature), ]
})

# Get LFC for genes in each of the 4 subgroups
genexp.comp <- lapply(names(hic.anno.list), function(gr) {
    splitcomp <- split(as.data.frame(gr), gr$transition)
})
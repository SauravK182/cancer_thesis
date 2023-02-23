#' @title Creates samplesheet for `DiffBind`/`ChIPQC`
#'
#' @description
#' Given a minimal `data.frame` object, and directories to `.bam` and `.xls` peak files,
#' constructs the sample sheet object used for the `DiffBind` and `ChIPQC` Bioconductor packages.
#' Note that as of the current development cycle, this function assumes peak calling was performed with
#' MACS, and that for ChIP-seq experiments, corresponding input controls were used. Please note that
#' several limiting assumptions are made based on the naming convention of your .bam and .xls files.
#' This function will automatically match bamReads and bamControls to the appropriate experiment, and
#' thus requires that naming conventions are uniform for all experiments, and that the user-supplied
#' information in `df` can be found in filenames (when requested - see `df` below)
#' @param df Minimal `data.frame` object. Assumes that the following columns are present:
#' * `SampleID` - a unique ID given to each experiment
#' * `Tissue` - Cell lines or biological samples used in the experiment. These strings should be present in your
#' `.bam` and `.xls` filenames.
#' * `Factor` - the protein/histone mark being ChIPped. This string should be in your `.bam` and `.xls` filenames.
#' For ATAC-seq experiments, you may use something like `atac-seq`, but this string should also be in the filenames.
#' * `Condition` - what condition the given experiment is in (e.g., WT, KO, primary tumor, metastasis, etc.)
#' * `Replicate` - replicate number for the experiment. This number should be present in your read/peak filenames.
#' @param bam.dir Directory where all `.bam` files can be found (should include ALL alignment files!)
#' @param peaks.dir Directory where all `.xls` peak files from MACS2 can be found (should include ALL peak files!)
#' @param atac Boolean indicating whether or not the experiments were ATAC-seq (i.e., have no input controls).
#' Default = FALSE
#'
#' @return Completed sample sheet that can be immediately used in `dba()` from `DiffBind` or alternative functions
#' in `ChIPQC`.
#' @export
#'
#' @examples
#' # Set up data frame for create.samplesheet()
#' tissue.cai <- rep(c("MDA-MB-231", "MDA-MB-231-BrM2", "MDA-MB-231-LM2"), each = 2)
#' factor.cai <- rep("h3k27ac", times = length(tissue.cai))
#' condition.cai <- rep(c("primary", "metastasis", "metastasis"), each = 2)
#' replicate.cai <- rep(c(1, 2), times = (length(tissue.cai) / 2))
#' sample.id.cai <- paste0(tissue.cai, replicate.cai)
#' cai.df <- data.frame(SampleID = sample.id.cai,
#'                   Tissue = tissue.cai,
#'                   Factor = factor.cai,
#'                   Condition = condition.cai,
#'                   Treatment = treatment.cai,
#'                   Replicate = replicate.cai)
#'
#' # Set up directories
#' cai.bam.dir <- file.path(getwd(), cai.data.dir, "chip-seq/aligned")
#' cai.peaks.dir <- file.path(getwd(), cai.data.dir, "chip-seq/called_peaks")

#' # Call samplesheet to use for dba()
#' cai.samplesheet <- create_samplesheet(cai.df, cai.bam.dir, cai.peaks.dir)
#' cai.dba <- dba(sampleSheet = cai.samplesheet)
create_samplesheet <- function(df, bam.dir, peaks.dir, atac = FALSE) {
    require(tidyverse)
    # Get .bam reads and .xls peak files based on the tissue, factor, and replicate
    bam.files <- c()
    peak.files <- c()
    for (i in seq_len(nrow(df))) {
        tissue <- df$Tissue[i]
        factor <- df$Factor[i]
        replicate <- df$Replicate[i]
        base.pattern <- paste0("(", tissue, "-|", tissue, "_)", factor, ".*", replicate, ".*")

        # Use base pattern to get each bam file for the tissue, factor, and replicate
        bam.name <- list.files(path = bam.dir, pattern = paste0(base.pattern, "bam$"))
        bam.file <- file.path(bam.dir, bam.name)
        bam.files <- c(bam.files, bam.file)

        # Use same approach as above to get the associated peak file
        peaks.name <- list.files(path = peaks.dir, pattern = paste0(base.pattern, "xls$"))
        peaks.file <- file.path(peaks.dir, peaks.name)
        peak.files <- c(peak.files, peaks.file)
    }
    if (! atac) {   # i.e., if there are control files
        # Search for files in bam.dir to isolate corresponding control treatments and add to sample sheet
        control.files <- list.files(path = bam.dir, pattern = ".*input.*bam$")  # gets list of all input control files
        # Order files with longest name first so that files w/ same basename are not counted twice (see below)
        control.files <- control.files[order(nchar(control.files), decreasing = TRUE)]
        control.treat.pairs <- list()
        control.id <- c()
        for (file in control.files) {
            treat.list <- c()
            basename <- gsub(pattern = "(.*\\/)", x = file, replacement = "") %>%
                        gsub(pattern = "[.]bam", replacement = "")      # removes leading directory and .bam, like basename in Bash
            # () denotes capturing group, treats multiple characters as one unit
            culture <- gsub(pattern = "-input.*|_input.*", replacement = "", x = basename)
            num.control <- gsub(pattern = ".*input", replacement = "", x = basename)
            for (trt.file in bam.files) {
                all.files <- c()
                for (vec in control.treat.pairs) {
                    all.files <- c(all.files, vec)  # to get a vector of all files already saved
                }
                if (trt.file %in% all.files) {    # if treat file is already in the list, skip

                } else if (grepl(pattern = paste0(culture, ".*", num.control, ".*"), x = trt.file)) {
                    treat.list <- c(treat.list, trt.file)   # if control corresponds to treatment, add
                }
            }
            if (length(treat.list) != 0) {
                control.treat.pairs[[file.path(bam.dir, file)]] <- treat.list
                control.id <- c(control.id, paste0(culture, num.control, "c"))
            }
        }
        # Now create the bamControl vector to go into the final spreadsheet
        bam.control <- c()
        control.id.reordered <- c()
        # Match the indices of the bam reads with the corresponding controls
        for (i in 1:length(control.treat.pairs)) {
            for (entry in control.treat.pairs[[i]]) {
                control.idx <- match(entry, bam.files)
                bam.control[control.idx] <- names(control.treat.pairs[i])
                control.id.reordered[control.idx] <- control.id[i]
            }
        }
        # Assemble final sheet
        sample.sheet <- cbind(df,
                            bamReads = bam.files,
                            ControlID = control.id.reordered,
                            bamControl = bam.control,
                            Peaks = peak.files,
                            PeakCaller = rep("macs", times = nrow(df)))
        return(sample.sheet)
    } else {
        sample.sheet <- cbind(df,
                              bamReads = bam.files,
                              Peaks = peak.files,
                              PeakCaller = rep("macs", times = nrow(df)))
        return(sample.sheet)
    }
}

#' @title Automate differential binding analysis with `DiffBind`
#'
#' @description
#' Takes in a samplesheet (or minimal `df` to be created into a samplsheet), list of contrasts, and a value for `summits`
#' along with other optional parameters to automate the repeated implementation of the `DiffBind` pipeline for differential
#' binding analysis
#' @param sample.df Either a minimal `data.frame` to be passed to `create_samplesheet()` or a sample sheet for use with `dba`.
#' If the `bamReads` column is present in `sample.df`, the function assumes the user has already set up the full sample sheet.
#' @param contrast.var Variable of interest (column in `sample.df`) for which DB analysis should be done as a main effect of.
#' @param summit.val Integer value to pass for `summits` in the `dba.count()` function. The `summits` parameter in `dba.count()` will
#' be used to set all consensus peaks to a fixed width after counting read pileup and recentering to the point of greatest
#' read density. The goal of `summits`, as described by `DiffBind` author Rory Stark, is to take a subsample of the merged peak
#' that is highly enriched for signal to be used as representative of the entire peak region. This mitigates the effect of background
#' noise (reads corresponding to nonspecific IP) on the construction of the affinity matrix.
#' @param formula.vec Character vector describing the formula design when `DiffBind` is passing the affinity matrix to `DESeq2`.
#' Default: equal to `contrast.var` (i.e., assumes the user is only interested in the effect of `contrast.var` on binding affinity).
#' @param contrasts List of character vectors describing all contrasts desired to be run by the user. Each character vector
#' should be of length 2 and correspond to entries in the column for `contrast.var`. Please note, the first-listed contrast
#' will be in the numerator of the LFC calculation.
#' Default: will run contrast on all possible unique combinations of the values in `contrast.var`
#' @param blacklist `DiffBind` blacklist object indicating what blacklist to enforce. Set to `FALSE` if user does not wish to
#' enforce any blacklists.
#' Default: DBA_BLACKLIST_GRCH37
#' @param atac Boolean value indicating whether the experiment corresponds to ChIP-seq data (`FALSE`) or ATAC-seq data (`TRUE`).
#' If the experiment is ATAC-seq data, the option `bUseSummarizeOverlaps` in `dba.count` will be set to `FALSE` to prevent
#' memory error crashes.
#' @param threshold FDR-adjusted p-value at which to call significantly differentially enriched regions.
#' Default: 0.05
#' @param lfc Real number defining the base LFC at which significantly differentially enriched regions are called.
#' Default: 0 (i.e., no LFC threshold)
#'
#' @return List of 3 objects:
#' * `dba.obj`, representing the original DBA object (with any enforced blacklists/greylists)
#' * `dba.counts`, representing the DBA object containing the count matrix for identified consensus peaks
#' * `dba.results`, list object containing the output of the `dba.analyze` function for each desired contrast.
#' @export
#'
#' @examples
#' test.analysis <- db_analysis(cai.samplesheet,
#'                              contrast.var = "Tissue",
#'                              summit.val = 500,
#'                              contrasts = list(c("MDA-MB-231-LM2", "MDA-MB-231")))
#' test.results.list <- test.analysis[[3]]
#' test.results <- test.results.list[[1]]
#' test.DB <- dba.report(test.results)
db_analysis <- function(sample.df,
                        contrast.var,
                        summit.val,
                        formula.vec = contrast.var,
                        contrasts = sample.df[, colnames(contrast.var) == contrast.var] %>% unique(),
                        blacklist = DBA_BLACKLIST_GRCH37,
                        atac = FALSE,
                        threshold = 0.05,
                        lfc = 0) {
    # Check if user supplied full samplesheet or base for sheet
    if ("data.frame" %in% class(sample.df)) {
        if ("bamReads" %in% colnames(sample.df)) {
            sample.sheet <- sample.df
        } else {
            sample.sheet <- create.samplesheet(sample.df)
        }
    } else {
        stop("Please supply a valid sample sheet data frame, or a data frame to be extended to a full
              sample sheet with create.samplesheet()")
    }

    # Get DBA object from sample sheet and apply blacklists and/or greylists
    dba.obj <- dba(sampleSheet = sample.sheet)
    if (! atac) {    # if not atac, apply greylist
        dba.obj <- dba.blacklist(dba.obj, blacklist = FALSE, greylist = TRUE)
    }
    if (blacklist != FALSE) {   # if blacklist not false, apply listed blacklist
        dba.obj <- dba.blacklist(dba.obj, blacklist = blacklist, greylist = FALSE)
    }

    # Get and save summary of peak widths after applying any blacklists
    peak.widths.sum <- summary(dba.obj$binding[, 3] - dba.obj$binding[, 2])

    # Count peaks using passed value for summits and normalize
    message(paste("Beginning counts for peak widths at a peak size of:", 2 * summit.val + 1, "bp"))
    if (2 * summit.val + 1 > peak.widths.sum[[2]]) {
        message <- paste("Selected summit size of", summit.val, "bp, resulting in a peak size of", 2 * summit.val + 1)
        message <- paste(message, "bp, is greater than the 1st quartile of total peak widths", paste0("(", peak.widths.sum[[2]], ")", "\n"))
        message <- paste0(message, "Summits will still be calculated as desired by user input.\n")
        message <- paste0(message, "Peak width may effect downstream analysis due to presence of background bases.")
        warning(message)
    }
    if (! atac) {
        dba.counts <- dba.count(dba.obj,
                                score = DBA_SCORE_RPKM,
                                summits = summit.val)
    } else {
        message("Now starting ATAC mode counting (bUseSummarizeOverlaps = FALSE)")
        dba.counts <- dba.count(dba.obj,
                                score = DBA_SCORE_RPKM,
                                summits = summit.val,
                                bUseSummarizeOverlaps = FALSE)  # to avoid memory errors, must set to FALSE
    }
    dba.norm <- dba.normalize(dba.counts,
                              method = DBA_DESEQ2,
                              normalize = DBA_NORM_NATIVE)

    # Get list of all possible contrasts - if user passed a list of contrasts, use these
    if (class(contrasts) == "list") {
        contrast.pairs <- contrasts
    } else {
        factorpairs.mat <- combn(contrasts, 2)
        contrast.pairs <- lapply(seq_len(ncol(factorpairs.mat)), function(x) factorpairs.mat[, x])
    }

    # Set up design formula from passed parameters
    collapsed.vec <- paste(formula.vec, collapse = "+") # will collapse vector "a b c" into "a+b+c"
    des.formula <- as.formula(paste0("~", collapsed.vec))

    # Set up contrasts and analyze
    dba.con.list <- lapply(1:length(contrast.pairs),
                           function(x) dba.contrast(dba.norm,
                                                    design = des.formula,
                                                    contrast = c(contrast.var, contrast.pairs[[x]])))
    dba.results <- lapply(dba.con.list, dba.analyze)

    # Return the DBA object, counts object, and the results list
    return(list(dba.obj, dba.counts, dba.results))
}


#' @title Simple wrapper for `annotatePeakInBatch`
#'
#' @description
#' Simple wrapper function for annotatePeakInBatch; will automatically take in a `GRanges` object and by default, annotate it
#' for the gene with the nearest TSS on both strands, assuming the peak falls within +/- 2kb of the TSS.
#'
#' @param dba.res Any `GRanges` object, with gene data, including for instance, a results object from `DiffBind`
#' @param annodata TSS or other gene data to use for `annotatePeakInBatch`. Default = `TSS.human.GRCh37`
#' @param output Specified mechanism for searching for nearest TSS from interval data. See `annotatePeakInBatch` for more.
#' Default = `nearestBiDirectionalPromoter`, which will report TSS's in both directions if the interval overlaps multiple TSS's.
#' @param region Annotation range for interval annotation. Default = `c(-2000, 2000)` - i.e., 2kb window around the peak
#' @param ... Any other parameter to be passed to `annotatePeakInBatch`
#'
#' @return Annotated GRanges object containing nearest TSS metadata for peaks within 2kb of one
#' @export
#'
#' @examples
anno_peak_gr37 <- function(dba.res, annodata = TSS.human.GRCh37,
                           output = "nearestBiDirectionalPromoters",
                           region = c(-2000, 2000), ...) {

    anno.peak <- annotatePeakInBatch(dba.res,
                                     AnnotationData = annodata,
                                     output = output,
                                     bindingRegion = region,
                                     ...)
    return(anno.peak)
}

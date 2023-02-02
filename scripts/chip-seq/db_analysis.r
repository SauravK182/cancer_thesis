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
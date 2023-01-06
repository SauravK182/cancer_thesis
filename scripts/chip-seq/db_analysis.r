create.samplesheet <- function(df, bam.dir, peaks.dir) {
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
            if (trt.file %in% control.treat.pairs) {    # if treat file is already in the list, skip

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
}

# Set up columns of the sampleSheet object for DiffBind
tissue <- rep(c("HPNE", "PANC-1", "Capan-1"), each = 2)
sample.id <- paste0(tissue, rep(c(1, 2), times = 3))
factor <- rep("H3K27ac", times = length(sample.id))
condition <- rep(c("normal", "primary", "metastasis"), each = 2)
treatment <- rep(NA, times = length(sample.id))
replicate <- rep(c(1, 2), times = 3)

test.df <- data.frame(SampleID = sample.id,
                      Tissue = tissue,
                      Factor = factor,
                      Condition = condition,
                      Treatment = treatment,
                      Replicate = replicate)
bam.dir <- file.path(main.data.dir, ren.data.dir, "chip-seq/aligned")
peaks.dir <- file.path(main.data.dir, ren.data.dir, "chip-seq/called_peaks")
test.sheet <- create.samplesheet(test.df, bam.dir, peaks.dir)

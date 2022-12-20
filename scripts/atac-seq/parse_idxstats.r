#!/usr/bin/env Rscript

# Use file("stdin") to read from stdin
input <- read.table(file("stdin"))

# Get rows where first col is non-canonical chrom (contains GL)
input.index <- grep("GL.*", input$V1)
input.subset <- input[input.index, ]

# For idxstats, 1st col is contig, 2nd col is length, 3rd col is num reads
non.canonical.length <- sum(input.subset$V2)
non.canonical <- sum(input.subset$V3)

# Create row entry for new df
non.canonical.df <- data.frame(V1 = "Non-canonical",
                               V2 = non.canonical.length,
                               V3 = non.canonical,
                               V4 = 0)

# Create new df with only one non-nuclear/canonical contig
nuclear.df <- rbind(input[-input.index, ], non.canonical.df)
write.table(nuclear.df, "nuclear_df.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
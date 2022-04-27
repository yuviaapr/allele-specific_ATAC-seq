#!/usr/bin/env Rscript

# Script to generate plot showing fragment size distribution of filtered mapped fragments.
# Usage: Rscript plot_fragment_size_distribution.R mapped_pairs.bam prefix

# Save path and prefix of the input bam file.

args <- commandArgs(trailingOnly = T)

bamFile <- args[1]

bamName <- args[2]

# Load library.

library(ATACseqQC)

# Generate and save plot.

pdf(paste0(bamName, "_fragment_size_distribution.pdf"), width = 8, height = 6)
	fragSizeDist(bamFile, bamName)
dev.off()



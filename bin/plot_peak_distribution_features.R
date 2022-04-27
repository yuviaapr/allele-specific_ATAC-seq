#!/usr/bin/env Rscript

# Script to generate plots with the distribution of peaks relative to genomic features.
# Usage: Rscript plot_peak_distribution_features.R peaks.bed prefix

# Save file path and prefix.

args <- commandArgs(trailingOnly = T)

peaks <- args[1]

plotName <- args[2]

# Load libraries.

library(ggupset)
library(org.Mm.eg.db)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# Load mouse gene database.

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Annotate peaks using the default definition of promoter (+-3 kb).
# Priority in genomic annotation for simple pie charts: Promoter > 5’ UTR > 3’ UTR > Exon > Intron > Downstream gene end > Intergenic.

peakAnno <- annotatePeak(peaks, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Mm.eg.db")

# Generate and save plots.

pdf(paste0(plotName, "_peak_annotation_distribution.pdf"), width = 12, height = 6)
	plotAnnoPie(peakAnno)
	upsetplot(peakAnno)
	vennpie(peakAnno)
dev.off()

# Save information about the annotations.

write.table(as.data.frame(peakAnno), paste0(plotName, "_peak_annotation_info.txt"), quote = F, row.names = F, sep = "\t")



# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load libraries
library(GenomicFeatures)
library(ChIPseeker)
library(tidyverse)
library(viridisLite)

# Load Snakemake parameters
bed.files <- snakemake@input[["bed"]]
gtf <- snakemake@input[["gtf"]]

# Load annotation database
txdb <- makeTxDbFromGFF(gtf)

# Add sample names to bed files
# Check which bed files are provided
if (any(grepl("\\.bed$", bed.files)) == TRUE) {
  samples <- sub(".*\\/([^\\/]+)\\.filtered.bed", "\\1", bed.files)
} else if (any(grepl("\\.narrowPeak$", bed.files)) == TRUE) {
  samples <- sub(".*\\/([^\\/]+)\\_peaks.narrowPeak", "\\1", bed.files)
} else if (any(grepl("\\.broadPeak$", bed.files)) == TRUE) {
  samples <- sub(".*\\/([^\\/]+)\\_peaks.broadPeak", "\\1", bed.files)
} 
names(bed.files) <- samples

# Annotate bed files
peakAnnoList <- lapply(bed.files,
                        annotatePeak,
                        TxDb = txdb,
                        tssRegion = c(-3000, 3000)
                        )

# Plot binding relative to TSS
pdf(snakemake@output[["dt"]],
    width = 10,
    height = length(bed.files) * 2.5)
plotDistToTSS(peakAnnoList,
              title =  "Distribution of binding sites relative to TSS") +
  theme(axis.line.y = element_line(linewidth = 0),
        axis.line.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = viridis(6),
                    name = "Distance to TSS") 
dev.off()

# Plot annotation bar
pdf(snakemake@output[["fd"]])
plotAnnoBar(peakAnnoList,
            title =  "Binding site distribution") +
  theme(axis.line.y = element_line(linewidth = 0),
        axis.line.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5))
dev.off()

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
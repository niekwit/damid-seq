# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load libraries
library(GenomicFeatures)
library(ChIPseeker)
library(tidyverse)

# Load Snakemake parameters
bed.files <- snakemake@input[["bed"]]
gtf <- snakemake@input[["gtf"]]

# Load annotation database
txdb <- makeTxDbFromGFF(gtf)

# Add sample names to bed files
samples <- sub(".*\\/([^\\/]+)\\.extended.bed", "\\1", bed.files)
names(bed.files) <- samples

# Annotate bed files
peakAnnoList <- lapply(bed.files,
                        annotatePeak,
                        TxDb = txdb,
                        tssRegion = c(-3000, 3000)
                        )

# Plot binding relative to TSS
pdf(snakemake@output[["tss"]])
plotDistToTSS(peakAnnoList,
              title =  "Distribution of peaks\nrelative to TSS") +
  theme(axis.line.y = element_line(size = 0),
        axis.line.x = element_line(size = 0.5),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 20))
dev.off()

# Plot annotation bar
pdf(snakemake@output[["bar"]])
plotAnnoBar(peakAnnoList) +
  theme(axis.line.y = element_line(size = 0),
        axis.line.x = element_line(size = 0.5),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        text = element_text(size = 20))
dev.off()

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
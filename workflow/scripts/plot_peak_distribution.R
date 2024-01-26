# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load libraries
library(ChIPpeakAnno)
library(GenomicFeatures)
library(stringr)

# Load Snakemake variables
bed.files <- snakemake@input[["bed"]]
gtf <- snakemake@input[["gtf"]]

# Create GRangesList object for all bed files
GRanges.list <- list()
for (b in bed.files) {
  name <- basename(str_replace(b, ".extended.bed", ""))
  GRanges.list[[name]] <- toGRanges(b,
                                    format = "BED",
                                    header = FALSE)
}
peaks <- GRangesList(GRanges.list)

# Load annotation data
txdb <- makeTxDbFromGFF(gtf)
annoData <- toGRanges(txdb, format = "gene")

# Save genomic element distribution object
pdf(snakemake@output[["ged"]])
genomicElementDistribution(peaks,
                           TxDb = txdb,
                           promoterRegion = c(upstream = 2000, 
                                              downstream = 500),
                           geneDownstream = c(upstream = 0, 
                                              downstream = 2000))
dev.off()

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
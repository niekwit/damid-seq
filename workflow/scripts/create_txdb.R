# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load libraries
library(GenomicFeatures)

# Load Snakemake parameters
gtf <- snakemake@input[["gtf"]]
txdb <- snakemake@output[["txdb"]]

# Create annotation database
txdb <- makeTxDbFromGFF(gtf)

save(txdb, file = txdb)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
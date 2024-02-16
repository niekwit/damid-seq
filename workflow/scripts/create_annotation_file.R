# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(rtracklayer)
library(tidyverse)

# Load Snakemake parameters
gtf <- snakemake@input[["gtf"]]

# Load GTF file
db <- rtracklayer::import(gtf)

# Extract relevant information
edb <- data.frame(geneId = db$gene_id, 
                  geneName = db$gene_name, 
                  geneBiotype = db$gene_biotype) %>%
  distinct()

# Save df to file as R object
save(edb, file = snakemake@output[["rdata"]])

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
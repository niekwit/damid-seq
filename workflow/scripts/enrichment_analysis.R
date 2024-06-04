# Redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

# Load libraries
library(tidyverse)
library(enrichR)
library(openxlsx)

# Load Snakemake variables
dbs <- snakemake@params[["dbs"]]
genome <- snakemake@params[["genome"]]
txt <- snakemake@input[["txt"]]
xlsx <- snakemake@output[["xlsx"]]

# Set up enrichr site
if (str_detect(genome, "dm")) {
  setEnrichrSite("FlyEnrichr")
} else if (str_detect(genome, "hg") || str_detect(genome, "mm")) {
  setEnrichrSite("Enrichr")
} else {
  stop("Genome not yet supported...")
}

# Load available databases
available_dbs <- listEnrichrDbs()

# Check if the requested databases are available
dbs_check <- dbs %in% available_dbs$libraryName
if (any(!dbs_check)) {
  print("Following databases are not available:")
  print(dbs[!dbs_check])
  print(paste("Available databases are:", paste(available_dbs$libraryName, collapse = ", ")))
  stop()
} else {
  print("All requested databases are available.")
}

# Get bound genes
genes <- read.csv(txt, header = FALSE, sep = " ") %>%
  pull(V2)

# Perform enrichment analysis
enrichment <- enrichr(genes, dbs)

# Write to excel file
write.xlsx(enrichment, xlsx)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
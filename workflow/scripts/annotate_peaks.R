# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load libraries
library(GenomicFeatures)
library(ChIPseeker)
library(tidyverse)

# Load Snakemake parameters
DB <- snakemake@input[["adb"]]
gtf <- snakemake@input[["gtf"]]
txt <- snakemake@output[["txt"]]
bed.file <- snakemake@input[["bed"]]

# Load annotation database
txdb <- makeTxDbFromGFF(gtf)

# Annotate bed file
peakAnno <- annotatePeak(bed.file,
                         tssRegion = c(-3000, 3000),
                         TxDb = txdb
                        )

# Tidy up annotation data
df <- as.data.frame(peakAnno@anno@elementMetadata@listData)
names(df)[1:3] <- c("peak_id", "fold_enrichment", "strand")

# Add gene names and gene biotype to annotation
load(DB)
df <- df %>%
  left_join(edb, by = "geneId")

# Write annotation to file
write.table(df,
            file = txt,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
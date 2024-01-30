# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load libraries
library(GenomicFeatures)
library(ChIPseeker)

# Load Snakemake parameters
bed <- snakemake@input[["bed"]]
gtf <- snakemake@input[["gtf"]]
txt <- snakemake@output[["txt"]]

# Load annotation database
txdb <- makeTxDbFromGFF(gtf)

# Annotate bed file
peakAnno <- annotatePeak(bed,
                         tssRegion = c(-3000, 3000),
                         TxDb = txdb)

# Tidy up annotation data
df <- as.data.frame(peakAnno@anno@elementMetadata@listData)
names(df)[1] <- "peak_id"
df <- df[, -c(2:4)]

# Write annotation to file
write.table(df,
            file = txt,
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
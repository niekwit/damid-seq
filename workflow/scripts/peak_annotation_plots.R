# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load libraries
library(GenomicFeatures)
library(ChIPseeker)

# Load Snakemake parameters
bed.files <- snakemake@input[["bed"]]
txdb <- snakemake@input[["txdb"]]

# Load annotation database
load(txdb)

# Annotation list of all bed files
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

# Feature distribution
pdf(snakemake@output[["fd"]])
plotAnnoBar(peakAnnoList)
dev.off()

# Distance to TSS
pdf(snakemake@output[["dt"]])
plotDistToTSS(peakAnnoList)
dev.off()

# Pathway analysis
genes <- lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) <- sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster = genes,
                           fun = "enrichKEGG",
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH")

pdf(snakemake@output[["pa"]])
dotplot(compKEGG,
        showCategory = 10,
        title = "KEGG Pathway Enrichment Analysis")
dev.off()

# Overlap of peaks and annotated genes
genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
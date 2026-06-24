# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(damidBind)
library(tidyverse)

# Load Snakemake parameters
bedgraphs <- snakemake@input[["bg"]]
outdir_bg <- snakemake@params[["outdir_bg"]]
peaks_gff <- snakemake@input[["gff"]]
outdir_peaks <- snakemake@params[["outdir_peaks"]]
norm_method <- snakemake@params[["norm_method"]]
genome <- snakemake@params[["genome"]]
fdr <- snakemake@params[["fdr"]]
filter_occupancy <- snakemake@params[["filter_occupancy"]]
organism <- case_when(
  str_detect(genome, "hg") ~ "homo sapiens",
  str_detect(genome, "mm") ~ "mus musculus",
  str_detect(genome, "dm") ~ "drosophila melanogaster",
  TRUE ~ "homo sapiens"
)
diff_diagn_plot <- snakemake@output[["diff_diagn_plot"]]

# Create output directories and symlink input files into them
# This is necessary because damidBind expects all files to be
# in the same directory
dir.create(outdir_bg, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_peaks, recursive = TRUE, showWarnings = FALSE)

print("Creating simlinks to original data in directory struture for damidBind")
for (bg in bedgraphs) {
  dir_name <- basename(dirname(bg))
  new_name <- sub(
    "^(.*?)(-vs-Dam.*)",
    paste0("\\1_", dir_name, "\\2"),
    basename(bg)
  )
  dest <- file.path(outdir_bg, new_name)
  if (!file.exists(dest)) file.symlink(normalizePath(bg), dest)
}

for (gff in peaks_gff) {
  dir_name <- basename(dirname(gff))
  new_name <- sub(
    "\\.peaks\\.gff$",
    paste0("_", dir_name, ".peaks.gff"),
    basename(gff)
  )
  dest <- file.path(outdir_peaks, new_name)
  if (!file.exists(dest)) file.symlink(normalizePath(gff), dest)
}

# Load data
print("Loading peak data and bedGraphs")
data <- load_data_peaks(
  binding_profiles_path = outdir_bg,
  peaks_path = outdir_peaks,
  norm_method = norm_method,
  plot_diagnostics = FALSE,
  organism = organism
)

# Prepare conditions for differential binding analysis
conditions <- c()
for (bg in bedgraphs) {
  condition <- sub("-vs-Dam.*", "", basename(bg))
  conditions <- c(conditions, condition)
}
conditions <- unique(conditions)

# Convert to named vector for damidBind
names(conditions) <- conditions

# Find differentially bound peaks
print("Finding differentially bound peaks")
pdf(diff_diagn_plot)
results <- tryCatch(
  differential_binding(
    data,
    cond = conditions,
    plot_diagnostics = TRUE,
    fdr = fdr,
    filter_occupancy = filter_occupancy
  ),
  error = function(e) {
    if (grepl("subscript out of bounds", conditionMessage(e))) {
      stop(
        "differential_binding() failed: too few loci passed the ",
        "occupancy filter. Try setting 'filter_occupancy' to an ",
        "integer lower than the number of replicates (e.g. 1) in ",
        "config.yaml under differential_peaks.",
        call. = FALSE
      )
    }
    stop(e)
  }
)
dev.off()

# Create Venn diagram of differentially bound peaks
pdf(snakemake@output[["venn"]])
plot_venn(results)
dev.off()

# Create volcano plot of differentially bound peaks
pdf(snakemake@output[["volcano"]])
plot_volcano(results, label_config = list(clean_names = TRUE))
dev.off()

# Save regions to CSV
print("Saving differentially bound peaks to CSV")
df <- analysisTable(results) %>%
  rownames_to_column("peak_id") %>%
  # split peak_id into chrom, start, end
  # keep peak_id for later use
  separate(
    peak_id,
    into = c("chrom", "start", "end"),
    sep = "[:-]",
    convert = TRUE,
    remove = FALSE
  )

write.csv(df, snakemake@output[["csv"]], row.names = FALSE)

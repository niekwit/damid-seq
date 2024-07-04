# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)

# Plot bar graph of FRiP values
total_read_counts <- snakemake@input[["total_read_counts"]]
peak_read_counts <- snakemake@input[["peak_read_count"]]

# df to store FRiP values
df <- data.frame(sample = character(),
                 read_in_peaks = numeric(),
                 total_reads = numeric(),
                 fraction_of_reads_in_peaks = numeric())

# Add FRiP values to df
for (i in seq_along(total_read_counts)) {
  sample <- snakemake@wildcards[["bg_sample"]][i]
  total_reads <- as.numeric(read_lines(total_read_counts[i]))
  peak_reads <- as.numeric(read_lines(peak_read_counts[i]))
  frip <- peak_reads / total_reads
  
  df <- df %>%
    add_row(sample = sample,
            read_in_peaks = peak_reads,
            total_reads = total_reads,
            fraction_of_reads_in_peaks = frip)
}

# Create plot
p <- ggplot(df, aes(x = sample, 
                    y = fraction_of_reads_in_peaks)) +
  geom_bar(stat = "identity", 
           fill = "#419179",
           colour = "black") +
  theme_cowplot(18) +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     limits = c(0, 1)) +
  labs(title = NULL,
       x = NULL,
       y = "Fraction of reads in peaks") +
  theme(plot.title = element_text(hjust = 0.5))

# Save plot
ggsave(snakemake@output[[1]], p)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
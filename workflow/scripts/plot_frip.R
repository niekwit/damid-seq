# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)

# Plot bar graph of FRiP values
total_read_counts <- snakemake@input[["total_read_count"]]
peak_read_counts <- snakemake@input[["peak_read_count"]]

# Create df to store FRiP values
df <- data.frame(sample = character(),
                 read_in_peaks = numeric(),
                 total_reads = numeric(),
                 fraction_of_reads_in_peaks = numeric())

# Add FRiP values to df
for (i in seq_along(total_read_counts)) {
  sample <- basename(total_read_counts[i]) %>% str_remove(".total.count")
  sample <- paste(basename(dirname(total_read_counts[i])), sample, sep="_")
  total_reads <- as.numeric(read_lines(total_read_counts[i]))
  peak_reads <- as.numeric(read_lines(peak_read_counts[i]))
  frip <- peak_reads / total_reads
  
  df <- df %>%
    add_row(sample = sample,
            read_in_peaks = peak_reads,
            total_reads = total_reads,
            fraction_of_reads_in_peaks = frip)
}

# Save df to file
write.csv(df, snakemake@output[["csv"]], row.names = FALSE)

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
       y = "Fraction of reads in peaks")

# Save plot
ggsave(snakemake@output[["plot"]], p)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
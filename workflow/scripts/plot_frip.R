# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)


# Read FRIP values from input file
csv <- snakemake@input[["csv"]]
df <- read.csv(csv, header = TRUE)

# Create plot
p <- ggplot(df, aes(x = name, y = fraction_in_peaks)) +
     geom_bar(stat = "identity", fill = "#419179", colour = "black") +
     theme_cowplot(18) +
     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
     scale_y_continuous(
          expand = expansion(mult = c(0, 0.1)),
          limits = c(0, 1)
     ) +
     labs(title = NULL, x = NULL, y = "Fraction of reads in peaks")

# Save plot
ggsave(snakemake@output[["plot"]], p)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")

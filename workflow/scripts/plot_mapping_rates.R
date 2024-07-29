# redirect R output to log
slog <- file(snakemake@log[[1]], open = "wt")
sink(slog, type = "output")
sink(slog, type = "message")

# Load libraries
library(tidyverse)
library(cowplot)

log.files <- snakemake@input[["log"]]
# Data frame to store mapping rates of all experiments
mapping.rates.all <- data.frame(sample = character(),
                                overall_mapping_rate = numeric())

# Extract mappings rates from each log file (one log file per experiment)
for (i in seq_along(log.files)) {
  # Get dir name (dir = experiment)
  dir <- basename(dirname(log.files[[i]]))

  # Get sample names and add dir name
  sample <- str_replace(paste0(dir, "_",basename(log.files[[i]])), ".log", "")

  # Read log file
  log <- readLines(log.files[[i]])
  
  # Extract mapping rate
  rate <- log[grepl("% overall alignment rate", log)]
  rate <- as.numeric(str_replace(rate, "% overall alignment rate", ""))

  # Extract mapping rates
  mapping.rates <- data.frame(sample = sample,
                              overall_mapping_rate = rate)

  # Add to data frame with all data
  mapping.rates.all <- rbind(mapping.rates.all, mapping.rates)
}

# Plot data and save
p <- ggplot(mapping.rates.all, 
            aes(x = sample, 
                y = overall_mapping_rate)) +
  geom_bar(stat = "identity",
           position = "dodge",
           colour = "black",
           fill = "#419179") +
  theme_cowplot(18) +
  theme(plot.margin = margin(t = 0.5,
                             r = 1.5,
                             b = 0.5,
                             l = 0.5,
                             unit = "cm")) +
  labs(x = NULL,
       y = "Overall alignment rates (%)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     limits = c(0, 100)) +
  scale_x_discrete(guide = guide_axis(angle = 45))

# Save to file
ggsave(snakemake@output[["pdf"]], p)

# Close redirection of output/messages
sink(slog, type = "output")
sink(slog, type = "message")
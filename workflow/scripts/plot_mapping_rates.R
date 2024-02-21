# redirect R output to log
slog <- file(snakemake@log[[1]], open = "wt")
sink(slog, type = "output")
sink(slog, type = "message")

# Load libraries
library(tidyverse)
library(cowplot)

# Get log files and dirs
# These log files are empty but use them to get 
# log dirs and find proper log files
log.files <- snakemake@input["log"]
log.dirs <- unlist(lapply(log.files, dirname))

# Get proper log files (log file starts with pipeline- and ends with .log)
log.files <- list.files(log.dirs,
                        pattern = "pipeline-.*.log",
                        full.names = TRUE)

# Data frame to store mapping rates of all experiments
mapping.rates.all <- data.frame(sample = character(),
                                overall_mapping_rate = numeric())

# Extract mappings rates from each log file (one log file per experiment)
for (i in seq_along(log.files)) {
  # Get dir name (dir = experiment)
  dir <- basename(log.dirs[[i]])

  # Read part of log file that contains mapping data
  log.section <- system(paste("sed '/Reading data files/,/Reading GATC file/!d'", log.files[[i]]), 
                        intern = TRUE)

  # Get sample names and add dir name
  sample.names <- log.section[grepl("Now working on ", log.section)]
  sample.names <- str_replace(sample.names, "Now working on ", "")
  sample.names <- paste0(dir, "_", str_replace(sample.names, " ...", ""))

  # Get line numbers where overall mapping rate is printed
  rate.lines <- grep("% overall alignment rate", log.section)
  
  # Extract mapping rates
  rates <- as.numeric(str_extract(log.section[rate.lines], "\\d+\\.\\d+"))

  # Extract mapping rates
  mapping.rates <- data.frame(sample = sample.names,
                              overall_mapping_rate = rates)

  # Add to data frame with all data
  mapping.rates.all <- rbind(mapping.rates.all, mapping.rates)
}

# Plot data and save
p <- ggplot(mapping.rates.all, aes(x = sample, y = overall_mapping_rate)) +
  geom_bar(stat = "identity",
           position = "dodge",
           colour = "black",
           fill = "aquamarine4") +
  theme_cowplot(18) +
  theme(plot.margin = margin(t = 0.5, r = 1.5, b = 0.5, l = 0.5, unit = "cm")) +
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

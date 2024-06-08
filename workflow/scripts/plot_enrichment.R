# Redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(openxlsx)
library(cowplot)

# Load Snakemake variables
xlsx <- snakemake@input[["xlsx"]]
dir_name <- snakemake@params[["dir_name"]]
terms <- snakemake@params[["terms"]]

# Create output directory
dir.create(dir_name, showWarnings = FALSE)

# Get sheet names
sheets <- getSheetNames(xlsx)

# Plot results for each sheet
for (i in sheets) {
  # Read data
  df <- read.xlsx(xlsx, sheet = i) %>%
    mutate(log.P.value = -log10(Adjusted.P.value)) %>%
    arrange(desc(log.P.value))
  
  # Convert character overlap ("n/m") to numerical ratio
  df$Ratio <- str_split(df$Overlap, "/") %>%
    map_dbl(~ as.numeric(.x[1]) / as.numeric(.x[2]))
  
  # Wrap terms longer than 60 characters over multiple lines
  df$Term <- str_wrap(df$Term, 60)
  
  # Relevel terms to avoid alphabetical sorting
  df$Term <- fct_rev(factor(df$Term, levels = df$Term))
  
  # Plot
  p <- ggplot(df[1:terms,], aes(x = log.P.value,
                                y = Term)) +
    geom_point(aes(fill = Ratio), 
               alpha = 1,
               size = 12,
               shape = 21,
               colour = "black") +
    scale_fill_gradient(low = "white", 
                        high = "forestgreen",
                        guide = guide_colorbar(frame.colour = "black", 
                                               ticks.colour = "black")) +
    theme_cowplot(18) +
    theme(axis.text.x = element_text(angle = 45, 
                                     hjust = 1)) +
    labs(x = "-log10(Adjusted P value)", 
         y = NULL) +
      ggtitle(i)
  
  # Save plot
  ggsave(paste0(dir_name, "/", i, ".pdf"), 
         p, 
         width = 12, 
         height = terms * 0.9)
}

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
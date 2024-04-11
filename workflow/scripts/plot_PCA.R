# redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# load required libraries
library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(ggrepel)
library(reshape2)
library(scales)

#### PCA plot ####
# Load PCA data
data <- read.delim(snakemake@input[[1]],
                   header = TRUE,
                   skip = 1)

# Unique sample conditions
colnames(data) <- gsub(".ext300", "", colnames(data))
samples <- colnames(data)[2:(ncol(data) - 1)]
samples <- unique(str_split(samples, "\\.", n = 2, simplify = TRUE)[, 2])

# Set colours for plotting
if (length(samples) == 1) {
  colours <- "#1B9E77"
} else if (length(samples) == 2) {
  colours <- c("#1B9E77", "#D95F02")
} else if (length(genotypes) > 2) {
  colours <- brewer.pal(length(samples), "Dark2")
}
names(colours) <- samples

# Keep only components 1 and 2, transpose and add sample information
df <- data[1:2, ] %>%
  dplyr::select(-c("Component", "Eigenvalue")) %>%
  t() %>%
  as.data.frame() %>%
  mutate(sample = rownames(.)) %>%
  rename(PC1 = 1,
         PC2 = 2) %>%
  mutate(sample_condition = str_split(sample, "\\.",
                                      n = 2,
                                      simplify = TRUE)[, 2],
         colour = colours[sample_condition])

# Calculate variance explained for each PC
PC1_var <- round((data$Eigenvalue[1] / sum(data$Eigenvalue)) * 100, 1)
PC2_var <- round((data$Eigenvalue[2] / sum(data$Eigenvalue)) * 100, 1)

# Create PCA plot
p <- ggplot(df,
            mapping = aes(x = PC1,
                          y = PC2,
                          colour = colour)) +
  geom_point(shape = 19,
             size = 5) +
  geom_label_repel(data = df,
                   aes(label = sample,
                       fill = NULL),
                   size = 5,
                   nudge_x = 0.5,
                   nudge_y = 0.5) +
  theme_cowplot(18) +
  labs(x = paste0("PC1: ", PC1_var, "% variance"),
       y = paste0("PC2: ", PC2_var, "% variance")) +
  theme(legend.position = "none")

# Save plot
ggsave(snakemake@output[["pca"]],
       p,
       height = 4,
       width = 6)

#### Scree plot ####
# Scale factor for utilising whole second y-axis range
# https://stackoverflow.com/questions/65559901/add-a-second-y-axis-to-ggplot
scalefactor <- max(data$Eigenvalue) / 100

# Prepare data for scree plot
df <- data %>%
  dplyr::select(c("Component", "Eigenvalue")) %>%
  mutate(Component = paste0("PC", Component)) %>%
  mutate(cumulative_variance = (cumsum(Eigenvalue) / sum(Eigenvalue) * 100 * scalefactor))

# Create scree plot
s <- ggplot(df, aes(Component, cumulative_variance)) +
  geom_bar(aes(Component, Eigenvalue),
           stat = "identity",
           colour = "black",
           fill = "aquamarine4") +
  geom_line(mapping = aes(x = Component,
                          y = cumulative_variance,
                          group = 1),
            colour = "red",
            linewidth = 1) +
  geom_point(mapping = aes(x = Component,
                           y = cumulative_variance),
             colour = "red",
             fill = "white",
             shape = 21,
             size = 6,
             stroke = 1.5) +
  theme_cowplot(16) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right = element_text(color = "red")) +
  scale_y_continuous(sec.axis = sec_axis(transform = ~ .x / scalefactor,
                                         breaks = seq(0, 100, 25),
                                         name = "Cumulative variance explained (%)"),
                     expand = expansion(mult = c(0, .05))) +
  labs(x = "Principal component",
       y = "Eigenvalue")

# Save plot
ggsave(snakemake@output[["scree"]],
       s,
       height = 4,
       width = 6)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
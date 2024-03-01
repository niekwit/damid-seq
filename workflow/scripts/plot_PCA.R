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

# Load sample information (replace any - with .)
sample_info <- read.csv("config/samples.csv", header = TRUE) %>%
  mutate(across(
    .cols = everything(),
    ~str_replace(., "-", ".")
  ))

# Get unique genotypes and treatments
genotypes <- unique(sample_info$genotype)
treatments <- unique(sample_info$treatment)

# Set colours (genotypes)
if (length(genotypes) < 3) {
  colours <- c("#1B9E77", "#D95F02")
} else {
  colours <- brewer.pal(length(genotypes), "Dark2")
}

# Set shapes (treatments)
shapes <- c(21, 22, 24, 23, 25)[seq_along(treatments)]

# Keep only components 1 and 2, transform and add sample information
df <- data[1:2, ] %>%
  dplyr::select(-c("Component", "Eigenvalue")) %>%
  t() %>%
  as.data.frame() %>%
  mutate(sample = rownames(.)) %>%
  rename(PC1 = 1,
         PC2 = 2) %>%
  left_join(sample_info, by = "sample")

# Calculate variance explained for each PC
PC1_var <- round((data$Eigenvalue[1] / sum(data$Eigenvalue)) * 100, 1)
PC2_var <- round((data$Eigenvalue[2] / sum(data$Eigenvalue)) * 100, 1)

# Create PCA plot
p <- ggplot(df,
            mapping = aes(x = PC1,
                          y = PC2,
                          fill = genotype,
                          shape = treatment)) +
  geom_point() +
  geom_label_repel(data = df,
                   aes(label = sample,
                       fill = NULL),
                   size = 5,
                   nudge_x = 0.5,
                   nudge_y = 0.5) +
  scale_fill_manual(values = colours) +
  scale_shape_manual(values = shapes) +
  theme_cowplot(16) +
  labs(x = paste0("PC1: ", PC1_var, "% variance"),
       y = paste0("PC2: ", PC2_var, "% variance"),
       Fill = "Genotype",
       shape = "Treatment")

# Save plot
ggsave(snakemake@output[["pca"]], p)


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
  scale_y_continuous(sec.axis = sec_axis(trans = ~ .x / scalefactor,
                                         breaks = seq(0, 100, 25),
                                         name = "Cumulative variance explained (%)"),
                     expand = expansion(mult = c(0, .05))) +
  labs(x = "Principal component",
       y = "Eigenvalue")

# Save plot
ggsave(snakemake@output[["scree"]], s)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")

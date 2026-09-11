# Species-level differential abundance analysis of Pseudomonas, Bradyrhizobium, Variovorax and Arthrobacter
# in the soybean rhizosphere at Leonard in 2022 using metagenomic data

# Load required packages
library(phyloseq)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(stringr)
library(openxlsx)
library(patchwork)

# Output directory
out_dir <- "Bacterial_species_level_differential_abundance_analysis"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}



# Pseudomonas

# Load bacterial metagenomics phyloseq object
ps <- readRDS("PS.metagenomics_bacteria.rds")

# Keep No fertilizer and High fertilizer treatments
ps_sub <- subset_samples(
  ps,
  Fertilizer %in% c("No fertilizer", "High fertilizer")
)

ps_sub <- prune_taxa(
  taxa_sums(ps_sub) > 0,
  ps_sub
)

sample_data(ps_sub)$Fertilizer <- factor(
  sample_data(ps_sub)$Fertilizer,
  levels = c("No fertilizer", "High fertilizer")
)

# Agglomerate taxa at the species level
ps_species <- tax_glom(
  ps_sub,
  taxrank = "Species",
  NArm = FALSE
)

ps_species <- prune_taxa(
  taxa_sums(ps_species) > 0,
  ps_species
)

# Taxonomy table
tx <- as.data.frame(tax_table(ps_species))
tx$feature <- rownames(tx)

# DESeq2 analysis
dds_species <- phyloseq_to_deseq2(
  ps_species,
  ~ Fertilizer
)

dds_species <- estimateSizeFactors(
  dds_species,
  type = "poscounts"
)

dds_species$Fertilizer <- relevel(
  dds_species$Fertilizer,
  ref = "No fertilizer"
)

dds_species <- DESeq(dds_species)

alpha <- 0.05

res_species <- results(
  dds_species,
  contrast = c(
    "Fertilizer",
    "High fertilizer",
    "No fertilizer"
  ),
  alpha = alpha
)

# Shrink log2 fold changes
res_species_shrunk <- lfcShrink(
  dds_species,
  contrast = c(
    "Fertilizer",
    "High fertilizer",
    "No fertilizer"
  ),
  res = res_species,
  type = "normal"
)

# Add taxonomy and keep Pseudomonas species
df_all <- as.data.frame(res_species_shrunk)
df_all$feature <- rownames(df_all)

df_pseudomonas <- df_all %>%
  left_join(
    tx[, c(
      "feature",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus",
      "Species"
    )],
    by = "feature"
  ) %>%
  filter(
    !is.na(Genus),
    str_trim(Genus) == "Pseudomonas",
    !is.na(log2FoldChange),
    !is.na(padj)
  ) %>%
  mutate(
    Species_clean = case_when(
      is.na(Species) | Species == "" ~ "sp.",
      grepl(
        "^(uncultured|unclassified|metagenome|bacterium)$",
        Species,
        ignore.case = TRUE
      ) ~ "sp.",
      TRUE ~ str_trim(as.character(Species))
    ),
    FullName = paste("Pseudomonas", Species_clean),
    Direction = ifelse(
      log2FoldChange >= 0,
      "Enriched in\nHigh fertilizer",
      "Enriched in\nNo fertilizer"
    ),
    Significant = padj < alpha
  )

# Keep significant Pseudomonas species and select up to 40 with largest absolute log2 fold changes
df_sig <- df_pseudomonas %>%
  filter(Significant == TRUE) %>%
  arrange(log2FoldChange) %>%
  slice_max(
    order_by = abs(log2FoldChange),
    n = 40
  )

if (nrow(df_sig) == 0) {
  stop("No significant Pseudomonas species found.")
}

# Save significant species results
write.xlsx(
  df_sig,
  file.path(
    out_dir,
    "DESeq2_high_vs_no_fertilizer_Pseudomonas_species_Leonard_2022_metagenomics.xlsx"
  ),
  overwrite = TRUE
)

# Calculate log10 baseMean
df_sig$log10baseMean <- log10(
  df_sig$baseMean + 1
)

# Set species order
df_sig$FullName <- factor(
  df_sig$FullName,
  levels = df_sig$FullName
)

# Log2 fold-change bar plot
p_bars <- ggplot(
  df_sig,
  aes(
    x = log2FoldChange,
    y = FullName,
    fill = Direction
  )
) +
  geom_bar(
    stat = "identity",
    color = "black",
    linewidth = 1.2,
    width = 1
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 1.2
  ) +
  scale_fill_manual(
    values = c(
      "Enriched in\nHigh fertilizer" = "red",
      "Enriched in\nNo fertilizer" = "darkgreen"
    )
  ) +
  scale_x_continuous(
    breaks = seq(-4, 4, 1)
  ) +
  labs(
    x = "log2FoldChange\n← No fertilizer          High fertilizer →",
    y = ""
  ) +
  theme_bw(base_size = 250) +
  theme(
    axis.line = element_line(color = "black", linewidth = 5),
    legend.position = "none",
    axis.text.y = element_text(size = 190, color = "black", face = "bold.italic"),
    axis.text.x = element_text(size = 250, color = "black", face = "bold"),
    axis.title.x = element_text(size = 250, face = "bold"),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 5)
  )

# BaseMean bar plot
p_basemean <- ggplot(
  df_sig,
  aes(
    x = log10baseMean,
    y = FullName
  )
) +
  geom_bar(
    stat = "identity",
    fill = "#0097A7",
    color = "black",
    linewidth = 2,
    width = 1
  ) +
  labs(
    x = "log10(baseMean + 1)",
    y = "",
    title = "Overall Abundance"
  ) +
  theme_bw(base_size = 250) +
  theme(
    axis.line = element_line(color = "black", linewidth = 5),
    legend.position = "none",
    plot.title = element_text(size = 250, face = "bold", hjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 250, face = "bold"),
    axis.title.x = element_text(size = 250, color = "black", face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 5),
    plot.margin = margin(10, 10, 10, 10)
  )

# Combine plots
final_species <- p_bars + p_basemean +
  plot_layout(
    widths = c(2.5, 1.0)
  )

# Dynamic figure height
n_species <- nrow(df_sig)
height_per_species <- 4
plot_height <- n_species * height_per_species

# Save figure
ggsave(
  filename = file.path(
    out_dir,
    "DESeq2_high_vs_no_fertilizer_Pseudomonas_species_Leonard_2022_metagenomics.pdf"
  ),
  plot = final_species,
  height = plot_height,
  width = 140,
  units = "in",
  bg = "white",
  limitsize = FALSE
)






# Bradyrhizobium 

# Load bacterial metagenomics phyloseq object
ps <- readRDS("PS.metagenomics_bacteria.rds")

# Keep No fertilizer and High fertilizer treatments
ps_sub <- subset_samples(
  ps,
  Fertilizer %in% c("No fertilizer", "High fertilizer")
)

ps_sub <- prune_taxa(
  taxa_sums(ps_sub) > 0,
  ps_sub
)

sample_data(ps_sub)$Fertilizer <- factor(
  sample_data(ps_sub)$Fertilizer,
  levels = c("No fertilizer", "High fertilizer")
)

# Agglomerate taxa at the species level
ps_species <- tax_glom(
  ps_sub,
  taxrank = "Species",
  NArm = FALSE
)

ps_species <- prune_taxa(
  taxa_sums(ps_species) > 0,
  ps_species
)

# Taxonomy table
tx <- as.data.frame(tax_table(ps_species))
tx$feature <- rownames(tx)

# DESeq2 analysis

dds_species <- phyloseq_to_deseq2(
  ps_species,
  ~ Fertilizer
)

dds_species <- estimateSizeFactors(
  dds_species,
  type = "poscounts"
)

dds_species$Fertilizer <- relevel(
  dds_species$Fertilizer,
  ref = "No fertilizer"
)

dds_species <- DESeq(dds_species)

alpha <- 0.05

res_species <- results(
  dds_species,
  contrast = c(
    "Fertilizer",
    "High fertilizer",
    "No fertilizer"
  ),
  alpha = alpha
)

# Shrink log2 fold changes
res_species_shrunk <- lfcShrink(
  dds_species,
  contrast = c(
    "Fertilizer",
    "High fertilizer",
    "No fertilizer"
  ),
  res = res_species,
  type = "normal"
)

# Add taxonomy and keep Bradyrhizobium species
df_all <- as.data.frame(res_species_shrunk)
df_all$feature <- rownames(df_all)

df_bradyrhizobium <- df_all %>%
  left_join(
    tx[, c(
      "feature",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus",
      "Species"
    )],
    by = "feature"
  ) %>%
  filter(
    !is.na(Genus),
    str_trim(Genus) == "Bradyrhizobium",
    !is.na(log2FoldChange),
    !is.na(padj)
  ) %>%
  mutate(
    Species_clean = case_when(
      is.na(Species) | Species == "" ~ "sp.",
      grepl(
        "^(uncultured|unclassified|metagenome|bacterium)$",
        Species,
        ignore.case = TRUE
      ) ~ "sp.",
      TRUE ~ str_trim(as.character(Species))
    ),
    FullName = paste("Bradyrhizobium", Species_clean),
    Direction = ifelse(
      log2FoldChange >= 0,
      "Enriched in\nHigh fertilizer",
      "Enriched in\nNo fertilizer"
    ),
    Significant = padj < alpha
  )

# Keep significant Bradyrhizobium species and select up to 40 with largest absolute log2 fold changes
df_sig <- df_bradyrhizobium %>%
  filter(Significant == TRUE) %>%
  arrange(log2FoldChange) %>%
  slice_max(
    order_by = abs(log2FoldChange),
    n = 40
  )

if (nrow(df_sig) == 0) {
  stop("No significant Bradyrhizobium species found.")
}

# Save significant species results
write.xlsx(
  df_sig,
  file.path(
    out_dir,
    "DESeq2_high_vs_no_fertilizer_Bradyrhizobium_species_Leonard_2022_metagenomics.xlsx"
  ),
  overwrite = TRUE
)

# Calculate log10 baseMean
df_sig$log10baseMean <- log10(
  df_sig$baseMean + 1
)

# Set species order
df_sig$FullName <- factor(
  df_sig$FullName,
  levels = df_sig$FullName
)

# Log2 fold-change bar plot
p_bars <- ggplot(
  df_sig,
  aes(
    x = log2FoldChange,
    y = FullName,
    fill = Direction
  )
) +
  geom_bar(
    stat = "identity",
    color = "black",
    linewidth = 1.2,
    width = 1
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 1.2
  ) +
  scale_fill_manual(
    values = c(
      "Enriched in\nHigh fertilizer" = "red",
      "Enriched in\nNo fertilizer" = "darkgreen"
    )
  ) +
  scale_x_continuous(
    breaks = seq(-4, 4, 1)
  ) +
  labs(
    x = "log2FoldChange\n← No fertilizer          High fertilizer →",
    y = ""
  ) +
  theme_bw(base_size = 250) +
  theme(
    axis.line = element_line(color = "black", linewidth = 5),
    legend.position = "none",
    axis.text.y = element_text(size = 190, color = "black", face = "bold.italic"),
    axis.text.x = element_text(size = 250, color = "black", face = "bold"),
    axis.title.x = element_text(size = 250, face = "bold"),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 5)
  )

# BaseMean bar plot
p_basemean <- ggplot(
  df_sig,
  aes(
    x = log10baseMean,
    y = FullName
  )
) +
  geom_bar(
    stat = "identity",
    fill = "#0097A7",
    color = "black",
    linewidth = 2,
    width = 1
  ) +
  labs(
    x = "log10(baseMean + 1)",
    y = "",
    title = "Overall Abundance"
  ) +
  theme_bw(base_size = 250) +
  theme(
    axis.line = element_line(color = "black", linewidth = 5),
    legend.position = "none",
    plot.title = element_text(size = 250, face = "bold", hjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 250, face = "bold"),
    axis.title.x = element_text(size = 250, color = "black", face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 5),
    plot.margin = margin(10, 10, 10, 10)
  )

# Combine plots
final_species <- p_bars + p_basemean +
  plot_layout(
    widths = c(2.5, 1.0)
  )

# Dynamic figure height

n_species <- nrow(df_sig)
height_per_species <- 4
plot_height <- n_species * height_per_species

# Save figure
ggsave(
  filename = file.path(
    out_dir,
    "DESeq2_high_vs_no_fertilizer_Bradyrhizobium_species_Leonard_2022_metagenomics.pdf"
  ),
  plot = final_species,
  height = plot_height,
  width = 140,
  units = "in",
  bg = "white",
  limitsize = FALSE
)




# Variovorax

# Load bacterial metagenomics phyloseq object
ps <- readRDS("PS.metagenomics_bacteria.rds")

# Keep No fertilizer and High fertilizer treatments
ps_sub <- subset_samples(
  ps,
  Fertilizer %in% c("No fertilizer", "High fertilizer")
)

ps_sub <- prune_taxa(
  taxa_sums(ps_sub) > 0,
  ps_sub
)

sample_data(ps_sub)$Fertilizer <- factor(
  sample_data(ps_sub)$Fertilizer,
  levels = c("No fertilizer", "High fertilizer")
)

# Agglomerate taxa at the species level
ps_species <- tax_glom(
  ps_sub,
  taxrank = "Species",
  NArm = FALSE
)

ps_species <- prune_taxa(
  taxa_sums(ps_species) > 0,
  ps_species
)

# Taxonomy table
tx <- as.data.frame(tax_table(ps_species))
tx$feature <- rownames(tx)

# DESeq2 analysis
dds_species <- phyloseq_to_deseq2(
  ps_species,
  ~ Fertilizer
)

dds_species <- estimateSizeFactors(
  dds_species,
  type = "poscounts"
)

dds_species$Fertilizer <- relevel(
  dds_species$Fertilizer,
  ref = "No fertilizer"
)

dds_species <- DESeq(dds_species)

alpha <- 0.05

res_species <- results(
  dds_species,
  contrast = c(
    "Fertilizer",
    "High fertilizer",
    "No fertilizer"
  ),
  alpha = alpha
)

# Shrink log2 fold changes
res_species_shrunk <- lfcShrink(
  dds_species,
  contrast = c(
    "Fertilizer",
    "High fertilizer",
    "No fertilizer"
  ),
  res = res_species,
  type = "normal"
)

# Add taxonomy and keep Variovorax species
df_all <- as.data.frame(res_species_shrunk)
df_all$feature <- rownames(df_all)

df_vario <- df_all %>%
  left_join(
    tx[, c(
      "feature",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus",
      "Species"
    )],
    by = "feature"
  ) %>%
  filter(
    !is.na(Genus),
    str_trim(Genus) == "Variovorax",
    !is.na(log2FoldChange),
    !is.na(padj)
  ) %>%
  mutate(
    Species_clean = case_when(
      is.na(Species) | Species == "" ~ "sp.",
      grepl(
        "^(uncultured|unclassified|metagenome|bacterium)$",
        Species,
        ignore.case = TRUE
      ) ~ "sp.",
      TRUE ~ str_trim(as.character(Species))
    ),
    FullName = paste("Variovorax", Species_clean),
    Direction = ifelse(
      log2FoldChange >= 0,
      "Enriched in\nHigh fertilizer",
      "Enriched in\nNo fertilizer"
    ),
    Significant = padj < alpha
  )

# Keep significant Variovorax species and select up to 40 with largest absolute log2 fold changes
df_sig <- df_vario %>%
  filter(Significant == TRUE) %>%
  arrange(log2FoldChange) %>%
  slice_max(
    order_by = abs(log2FoldChange),
    n = 40
  )

if (nrow(df_sig) == 0) {
  stop("No significant Variovorax species found.")
}

# Save significant species results
write.xlsx(
  df_sig,
  file.path(
    out_dir,
    "DESeq2_high_vs_no_fertilizer_Variovorax_species_Leonard_2022_metagenomics.xlsx"
  ),
  overwrite = TRUE
)

# Calculate log10 baseMean
df_sig$log10baseMean <- log10(
  df_sig$baseMean + 1
)

# Set species order

df_sig$FullName <- factor(
  df_sig$FullName,
  levels = df_sig$FullName
)

# Log2 fold-change bar plot
p_bars <- ggplot(
  df_sig,
  aes(
    x = log2FoldChange,
    y = FullName,
    fill = Direction
  )
) +
  geom_bar(
    stat = "identity",
    color = "black",
    linewidth = 1.2,
    width = 1
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 1.2
  ) +
  scale_fill_manual(
    values = c(
      "Enriched in\nHigh fertilizer" = "red",
      "Enriched in\nNo fertilizer" = "darkgreen"
    )
  ) +
  scale_x_continuous(
    breaks = seq(-4, 4, 1)
  ) +
  labs(
    x = "log2FoldChange\n← No fertilizer          High fertilizer →",
    y = ""
  ) +
  theme_bw(base_size = 250) +
  theme(
    axis.line = element_line(color = "black", linewidth = 5),
    legend.position = "none",
    axis.text.y = element_text(size = 190, color = "black", face = "bold.italic"),
    axis.text.x = element_text(size = 250, color = "black", face = "bold"),
    axis.title.x = element_text(size = 250, face = "bold"),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 5)
  )

# BaseMean bar plot
p_basemean <- ggplot(
  df_sig,
  aes(
    x = log10baseMean,
    y = FullName
  )
) +
  geom_bar(
    stat = "identity",
    fill = "#0097A7",
    color = "black",
    linewidth = 2,
    width = 1
  ) +
  labs(
    x = "log10(baseMean + 1)",
    y = "",
    title = "Overall Abundance"
  ) +
  theme_bw(base_size = 250) +
  theme(
    axis.line = element_line(color = "black", linewidth = 5),
    legend.position = "none",
    plot.title = element_text(size = 250, face = "bold", hjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 250, face = "bold"),
    axis.title.x = element_text(size = 250, color = "black", face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 5),
    plot.margin = margin(10, 10, 10, 10)
  )

# Combine plots
final_species <- p_bars + p_basemean +
  plot_layout(
    widths = c(2.5, 1.0)
  )

# Dynamic figure height
n_species <- nrow(df_sig)
height_per_species <- 4
plot_height <- n_species * height_per_species

# Save figure
ggsave(
  filename = file.path(
    out_dir,
    "DESeq2_high_vs_no_fertilizer_Variovorax_species_Leonard_2022_metagenomics.pdf"
  ),
  plot = final_species,
  height = plot_height,
  width = 140,
  units = "in",
  bg = "white",
  limitsize = FALSE
)







# Arthrobacter 

# Load bacterial metagenomics phyloseq object
ps <- readRDS("PS.metagenomics_bacteria.rds")

# Keep No fertilizer and High fertilizer treatments
ps_sub <- subset_samples(
  ps,
  Fertilizer %in% c("No fertilizer", "High fertilizer")
)

ps_sub <- prune_taxa(
  taxa_sums(ps_sub) > 0,
  ps_sub
)

sample_data(ps_sub)$Fertilizer <- factor(
  sample_data(ps_sub)$Fertilizer,
  levels = c("No fertilizer", "High fertilizer")
)

# Agglomerate taxa at the species level
ps_species <- tax_glom(
  ps_sub,
  taxrank = "Species",
  NArm = FALSE
)

ps_species <- prune_taxa(
  taxa_sums(ps_species) > 0,
  ps_species
)

# Taxonomy table
tx <- as.data.frame(tax_table(ps_species))
tx$feature <- rownames(tx)

# DESeq2 analysis
dds_species <- phyloseq_to_deseq2(
  ps_species,
  ~ Fertilizer
)

dds_species <- estimateSizeFactors(
  dds_species,
  type = "poscounts"
)

dds_species$Fertilizer <- relevel(
  dds_species$Fertilizer,
  ref = "No fertilizer"
)

dds_species <- DESeq(dds_species)

alpha <- 0.05

res_species <- results(
  dds_species,
  contrast = c(
    "Fertilizer",
    "High fertilizer",
    "No fertilizer"
  ),
  alpha = alpha
)

# Shrink log2 fold changes
res_species_shrunk <- lfcShrink(
  dds_species,
  contrast = c(
    "Fertilizer",
    "High fertilizer",
    "No fertilizer"
  ),
  res = res_species,
  type = "normal"
)

# Add taxonomy and keep Arthrobacter species
df_all <- as.data.frame(res_species_shrunk)
df_all$feature <- rownames(df_all)

df_arthrobacter <- df_all %>%
  left_join(
    tx[, c(
      "feature",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus",
      "Species"
    )],
    by = "feature"
  ) %>%
  filter(
    !is.na(Genus),
    str_trim(Genus) == "Arthrobacter",
    !is.na(log2FoldChange),
    !is.na(padj)
  ) %>%
  mutate(
    Species_clean = case_when(
      is.na(Species) | Species == "" ~ "sp.",
      grepl(
        "^(uncultured|unclassified|metagenome|bacterium)$",
        Species,
        ignore.case = TRUE
      ) ~ "sp.",
      TRUE ~ str_trim(as.character(Species))
    ),
    FullName = paste("Arthrobacter", Species_clean),
    Direction = ifelse(
      log2FoldChange >= 0,
      "Enriched in\nHigh fertilizer",
      "Enriched in\nNo fertilizer"
    ),
    Significant = padj < alpha
  )

# Keep significant Arthrobacter species and select up to 40 with largest absolute log2 fold changes
df_sig <- df_arthrobacter %>%
  filter(Significant == TRUE) %>%
  arrange(log2FoldChange) %>%
  slice_max(
    order_by = abs(log2FoldChange),
    n = 40
  )

if (nrow(df_sig) == 0) {
  stop("No significant Arthrobacter species found.")
}

# Save significant species results
write.xlsx(
  df_sig,
  file.path(
    out_dir,
    "DESeq2_high_vs_no_fertilizer_Arthrobacter_species_Leonard_2022_metagenomics.xlsx"
  ),
  overwrite = TRUE
)

# Calculate log10 baseMean
df_sig$log10baseMean <- log10(
  df_sig$baseMean + 1
)

# Set species order

df_sig$FullName <- factor(
  df_sig$FullName,
  levels = df_sig$FullName
)

# Log2 fold-change bar plot
p_bars <- ggplot(
  df_sig,
  aes(
    x = log2FoldChange,
    y = FullName,
    fill = Direction
  )
) +
  geom_bar(
    stat = "identity",
    color = "black",
    linewidth = 1.2,
    width = 1
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 1.2
  ) +
  scale_fill_manual(
    values = c(
      "Enriched in\nHigh fertilizer" = "red",
      "Enriched in\nNo fertilizer" = "darkgreen"
    )
  ) +
  scale_x_continuous(
    breaks = seq(-4, 4, 1)
  ) +
  labs(
    x = "log2FoldChange\n← No fertilizer          High fertilizer →",
    y = ""
  ) +
  theme_bw(base_size = 250) +
  theme(
    axis.line = element_line(color = "black", linewidth = 5),
    legend.position = "none",
    axis.text.y = element_text(size = 190, color = "black", face = "bold.italic"),
    axis.text.x = element_text(size = 250, color = "black", face = "bold"),
    axis.title.x = element_text(size = 250, face = "bold"),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 5)
  )

# BaseMean bar plot
p_basemean <- ggplot(
  df_sig,
  aes(
    x = log10baseMean,
    y = FullName
  )
) +
  geom_bar(
    stat = "identity",
    fill = "#0097A7",
    color = "black",
    linewidth = 2,
    width = 1
  ) +
  labs(
    x = "log10(baseMean + 1)",
    y = "",
    title = "Overall Abundance"
  ) +
  theme_bw(base_size = 250) +
  theme(
    axis.line = element_line(color = "black", linewidth = 5),
    legend.position = "none",
    plot.title = element_text(size = 250, face = "bold", hjust = 0.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 250, face = "bold"),
    axis.title.x = element_text(size = 250, color = "black", face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 5),
    plot.margin = margin(10, 10, 10, 10)
  )

# Combine plots
final_species <- p_bars + p_basemean +
  plot_layout(
    widths = c(2.5, 1.0)
  )

# Dynamic figure height
n_species <- nrow(df_sig)
height_per_species <- 4
plot_height <- n_species * height_per_species

# Save figure
ggsave(
  filename = file.path(
    out_dir,
    "DESeq2_high_vs_no_fertilizer_Arthrobacter_species_Leonard_2022_metagenomics.pdf"
  ),
  plot = final_species,
  height = plot_height,
  width = 140,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

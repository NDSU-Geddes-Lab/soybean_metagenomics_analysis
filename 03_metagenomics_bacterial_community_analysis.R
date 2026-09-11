# Soybean rhizosphere bacterial community analysis across different fertilizer treatments 
# at Leonard in 2022 using metagenomics data

# Load required packages
library(phyloseq)
library(vegan)
library(pairwiseAdonis)
library(ggplot2)
library(cowplot)
library(rstatix)
library(dplyr)
library(forcats)
library(writexl)
library(openxlsx)
library(ggh4x)
library(grid)
library(DESeq2)
library(data.tree)
library(ape)
library(ggtree)
library(patchwork)

# Output directory
out_dir <- "Metagenomics_bacterial_community_analysis_Leonard_2022"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Fertilizer order
fert_order <- c("No fertilizer", "Medium fertilizer", "High fertilizer")


# Bacterial Shannon diversity analysis across fertilizer treatments in rhizosphere at Leonard in 2022

# Load bacterial metagenomics phyloseq object
ps <- readRDS("PS.metagenomics_bacteria.rds")

f <- sample_data(ps)$Fertilizer |> as.character() |> trimws()
sample_data(ps)$Fertilizer <- fct_relevel(factor(f), fert_order)

# Calculate Shannon diversity
alpha_df <- estimate_richness(ps, measures = "Shannon")
alpha_df$SampleID <- rownames(alpha_df)
meta <- as(sample_data(ps), "data.frame")
meta$SampleID <- rownames(meta)
alpha_df <- left_join(alpha_df, meta, by = "SampleID")
alpha_df$Fertilizer <- fct_relevel(factor(trimws(as.character(alpha_df$Fertilizer))), fert_order)

# Shapiro-Wilk normality test
sw_results <- do.call(rbind, lapply(fert_order, function(g) {
  sub <- alpha_df[alpha_df$Fertilizer == g, ]
  sw <- shapiro.test(sub$Shannon)
  data.frame(
    Fertilizer = g,
    n = nrow(sub),
    W = round(sw$statistic, 4),
    p = round(sw$p.value, 4),
    Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
}))

print(sw_results)

# Normality not met: Kruskal-Wallis + pairwise Wilcoxon with BH correction
# Kruskal-Wallis test
kw_result <- alpha_df %>% kruskal_test(Shannon ~ Fertilizer)
print(as.data.frame(kw_result))

# Pairwise Wilcoxon test with BH correction
pwc_all <- alpha_df %>%
  wilcox_test(Shannon ~ Fertilizer, p.adjust.method = "BH", paired = FALSE) %>%
  add_significance("p.adj")

print(as.data.frame(pwc_all))

# Save statistics
kw_export <- as.data.frame(kw_result) %>%
  mutate(test = "Kruskal-Wallis", group1 = NA_character_, group2 = NA_character_, p.adj = NA_real_, p.adj.signif = NA_character_) %>%
  select(test, group1, group2, statistic, p, p.adj, p.adj.signif)

pwc_export <- as.data.frame(pwc_all) %>%
  mutate(test = "Pairwise Wilcoxon (BH adjusted)") %>%
  select(test, group1, group2, statistic, p, p.adj, p.adj.signif)

write_xlsx(
  bind_rows(kw_export, pwc_export),
  file.path(out_dir, "Alpha_diversity_statistics_rhizosphere_bacteria_Leonard_2022_metagenomics.xlsx")
)




# Non-metric Multidimensional Scaling or NMDS analysis using Bray-Curtis distance of rhizosphere 
# bacterial communities across fertilizer treatments at Leonard in 2022 

# Convert to relative abundance
set.seed(777)
ps <- readRDS("PS.metagenomics_bacteria.rds")
ps_prop_metag_b <- transform_sample_counts(ps, function(otu) otu / sum(otu))

# Bray-Curtis NMDS
ord.nmds.bray.metag.b <- ordinate(
  ps_prop_metag_b,
  method = "NMDS",
  distance = "bray"
)

# Fertilizer colors
custom_colors <- c(
  "No fertilizer" = "#8B1C62",
  "Medium fertilizer" = "#00008B",
  "High fertilizer" = "#009E73"
)

# Plot settings
POINT_SIZE <- 45
POINT_STROKE <- 3.5
JITTER_W <- 0.03
JITTER_H <- 0.03
TITLE_SIZE <- 130
AXIS_TITLE <- 130
AXIS_TEXT <- 140
FIG_W <- 38
FIG_H <- 25
FIG_DPI <- 600

# NMDS coordinates
ord_df <- plot_ordination(
  ps_prop_metag_b,
  ord.nmds.bray.metag.b,
  color = "Fertilizer",
  justDF = TRUE
)

# NMDS plot
nmds.plot <- ggplot(
  ord_df,
  aes(x = NMDS1, y = NMDS2, color = Fertilizer)
) +
  theme_cowplot() +
  geom_point(
    size = POINT_SIZE,
    stroke = POINT_STROKE,
    position = position_jitter(width = JITTER_W, height = JITTER_H)
  ) +
  scale_color_manual(values = custom_colors) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme(
    axis.line = element_line(linewidth = 2),
    legend.position = "right",
    plot.title = element_text(size = TITLE_SIZE, face = "bold", hjust = 0.5),
    axis.title = element_text(size = AXIS_TITLE, face = "bold"),
    axis.text = element_text(size = AXIS_TEXT, face = "bold")
  )

ggsave(
  filename = file.path(out_dir, "NMDS_rhizosphere_bacteria_Leonard_2022_metagenomics.png"),
  plot = nmds.plot,
  width = FIG_W,
  height = FIG_H,
  units = "in",
  dpi = FIG_DPI,
  bg = "white"
)

# PERMANOVA and pairwise PERMANOVA analysis of rhizosphere bacterial communities 
# across fertilizer treatments at Leonard in 2022 

bray_metag_b <- phyloseq::distance(ps_prop_metag_b, method = "bray")
meta_metag_b <- data.frame(sample_data(ps_prop_metag_b))

# Check dispersion across fertilizer treatments
set.seed(777)
disp_metag_b <- betadisper(bray_metag_b, meta_metag_b$Fertilizer)
dispersion_result <- anova(disp_metag_b)

print(dispersion_result)

# PERMANOVA across fertilizer treatments
set.seed(777)
permanova_result <- adonis2(
  bray_metag_b ~ Fertilizer,
  data = meta_metag_b,
  permutations = 999
)

print(permanova_result)

# Pairwise PERMANOVA across fertilizer treatments
set.seed(777)
pairwise_result <- pairwise.adonis2(
  bray_metag_b ~ Fertilizer,
  data = meta_metag_b,
  permutations = 999,
  p.adjust.m = "bonferroni"
)

print(pairwise_result)

# Save PERMANOVA and pairwise PERMANOVA results
pairwise_df <- do.call(rbind, lapply(names(pairwise_result), function(pair) {
  x <- pairwise_result[[pair]]
  if (!is.data.frame(x)) return(NULL)
  cbind(Comparison = pair, x)
}))

results <- list(
  Fertilizer_pairwise_PERMANOVA = pairwise_df,
  Fertilizer_PERMANOVA = permanova_result
)

write.xlsx(
  results,
  file.path(out_dir, "PERMANOVA_and_pairwise_PERMANOVA_rhizosphere_bacteria_Leonard_2022_metagenomics.xlsx"),
  overwrite = TRUE
)


# Stacked barplot showing relative abundance of rhizosphere bacterial genera (top 40)
# across fertilizer treatments at Leonard in 2022 

# Load bacterial metagenomics phyloseq object
ps <- readRDS("PS.metagenomics_bacteria.rds")

# Keep taxa with genus-level classification
ps.genus <- subset_taxa(ps, !is.na(Genus) & Genus != "")

# Agglomerate taxa at the genus level
ps.genus <- tax_glom(ps.genus, taxrank = "Genus", NArm = FALSE)

# Keep the 40 most abundant genera
top40 <- names(sort(taxa_sums(ps.genus), decreasing = TRUE)[1:40])
ps.genus <- prune_taxa(top40, ps.genus)

# Convert counts to relative abundance
ps.rel_amp <- transform_sample_counts(ps.genus, function(x) x / sum(x))

# Sample names
if (!"Sample.names" %in% colnames(sample_data(ps.rel_amp))) {
  sample_data(ps.rel_amp)$Sample.names <- sample_names(ps.rel_amp)
}

# Fertilizer facets
sample_data(ps.rel_amp)$FertilizerFacet <- factor(
  as.character(sample_data(ps.rel_amp)$Fertilizer),
  levels = c("No fertilizer", "Medium fertilizer", "High fertilizer")
)

# Color palette
my_colors <- c(
  "#C5CAE9", "#CDDC39", "#3949AB", "#1E88E5", "#00ACC1", "#00897B", "darkgreen", "#C0CA33", "#D7CCC8",
  "#FB8C00", "#F4511E", "#6D4C41", "#8D6E63", "#D84315", "#D81B60", "#5E35B1", "#0277BD", "#00838F",
  "#00B8D4", "#00BFA5", "#7CB342", "#757575", "#FFE0B2", "#FFB300", "#D1C4E9", "#757575", "#FDD835",
  "#AD1457", "#78909C", "#4527A0", "#283593", "#0277BD", "#00838F", "#00695C", "#2E7D32", "#9E9D24",
  "#F9A825", "#F8BBD0", "#546E7A", "#4E342E"
)

# Stacked relative abundance plot
p <- plot_bar(ps.rel_amp, x = "Sample.names", fill = "Genus") +
  ggh4x::facet_nested(
    ~ FertilizerFacet,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(values = my_colors) +
  labs(x = "Sample", y = "Relative abundance") +
  theme_bw(base_size = 26) +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 44, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_line(),
    axis.text.y = element_text(size = 44),
    axis.title = element_text(size = 44, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(0, "pt"),
    ggh4x.facet.nestline = element_line(linewidth = 2)
  )

# Save plot
ggsave(
  file.path(out_dir, "Stacked_barplot_top40_bacterial_genus_fertilizerlevels_Leonard_rhizosphere_2022_metagenomics.png"),
  p, width = 24, height = 12, dpi = 600, bg = "white")

# Save genus legend separately
p <- plot_bar(ps.rel_amp, x = "Sample.names", fill = "Genus") +
  scale_fill_manual(values = my_colors) +
  theme_bw(base_size = 26)

get_only_legend <- function(myplot) {
  tmp <- ggplot_gtable(ggplot_build(myplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

legend_only <- get_only_legend(p)

ggsave(file.path(out_dir, "Legend_top40_bacterial_genus_fertilizerlevels_Leonard_rhizosphere_2022_metagenomics.png"),
  plot = legend_only, width = 8, height = 12, units = "in", dpi = 600, bg = "white")





# Differential abundance analysis by DESeq2 of rhizosphere bacterial genera between
# the no fertilizer and high fertilizer treatments at Leonard in 2022

# Load bacterial metagenomics phyloseq object
ps <- readRDS("PS.metagenomics_bacteria.rds")

# Keep No and High fertilizer treatments
ps_sub <- subset_samples(ps, Fertilizer %in% c("No fertilizer", "High fertilizer"))
ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)

sample_data(ps_sub)$Fertilizer <- factor(
  sample_data(ps_sub)$Fertilizer,
  levels = c("No fertilizer", "High fertilizer")
)

# Agglomerate taxa at the genus level
ps_genus <- tax_glom(ps_sub, taxrank = "Genus", NArm = TRUE)

tx <- as.data.frame(tax_table(ps_genus))
ps_genus <- prune_taxa(!is.na(tx$Genus) & nzchar(tx$Genus), ps_genus)

tx <- as.data.frame(tax_table(ps_genus))
tx$feature <- rownames(tx)

# DESeq2 analysis
dds_genus <- phyloseq_to_deseq2(ps_genus, ~ Fertilizer)
dds_genus <- estimateSizeFactors(dds_genus, type = "poscounts")
dds_genus$Fertilizer <- relevel(dds_genus$Fertilizer, ref = "No fertilizer")
dds_genus <- DESeq(dds_genus)

alpha <- 0.05

res_genus <- results(
  dds_genus,
  contrast = c("Fertilizer", "High fertilizer", "No fertilizer"),
  alpha = alpha
)

res_genus_shrunk <- lfcShrink(
  dds_genus,
  contrast = c("Fertilizer", "High fertilizer", "No fertilizer"),
  res = res_genus,
  type = "normal"
)

# Keep significant genera
sig_g <- as.data.frame(res_genus_shrunk)
sig_g$feature <- rownames(sig_g)

sig_g <- sig_g %>%
  left_join(tx, by = "feature") %>%
  filter(
    !is.na(padj),
    padj < alpha,
    Kingdom == "Bacteria",
    !is.na(Phylum),
    Phylum != "",
    !is.na(Genus),
    Genus != "",
    Genus != "Unclassified"
  )

sig_g$Genus_clean <- make.unique(as.character(sig_g$Genus))

# Fill missing taxonomy for tree construction
for (col in c("Kingdom", "Phylum", "Class", "Order", "Family")) {
  sig_g[[col]][is.na(sig_g[[col]]) | sig_g[[col]] == ""] <- "Unclassified"
}

sig_g <- sig_g %>%
  arrange(Phylum, Class, Order, Family, Genus_clean)

# Build taxonomy tree
tax_df_g <- data.frame(
  Kingdom = sig_g$Kingdom,
  Phylum = sig_g$Phylum,
  Class = sig_g$Class,
  Order = sig_g$Order,
  Family = sig_g$Family,
  Genus_clean = sig_g$Genus_clean
)

tax_df_g$pathString <- paste(
  "Root",
  tax_df_g$Kingdom,
  tax_df_g$Phylum,
  tax_df_g$Class,
  tax_df_g$Order,
  tax_df_g$Family,
  tax_df_g$Genus_clean,
  sep = "/"
)

tax_df_g <- tax_df_g[!duplicated(tax_df_g$pathString), ]

tax_tree_g <- as.Node(tax_df_g)
writeLines(ToNewick(tax_tree_g), "temp_genus_tree.nwk")
phylo_tree_g <- read.tree("temp_genus_tree.nwk")

# Phylum colors
phylum_colors <- c(
  "Pseudomonadota" = "#4A148C",
  "Actinomycetota" = "#C2185B",
  "Bacillota" = "#1B5E20",
  "Bacteroidota" = "#B8860B",
  "Planctomycetota" = "#FDD835",
  "Chloroflexota" = "#FF6600",
  "Cyanobacteriota" = "#1565C0",
  "Chlamydiota" = "#78909C",
  "Thermotogota" = "#00CCCC",
  "Verrucomicrobiota" = "#F48FB1"
)

# Prepare tree-tip taxonomy
tip_tax_g <- data.frame(
  label = sig_g$Genus_clean,
  Phylum = sig_g$Phylum
)

tip_tax_g <- tip_tax_g[!duplicated(tip_tax_g$label), ]

valid_tips <- intersect(phylo_tree_g$tip.label, tip_tax_g$label)
phylo_tree_g <- keep.tip(phylo_tree_g, valid_tips)

# Get tree tip order
p_tree_tmp <- ggtree(phylo_tree_g) %<+% tip_tax_g

tree_order_g <- p_tree_tmp$data %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label)

sig_g <- sig_g[sig_g$Genus_clean %in% tree_order_g, ]
sig_g$Genus_clean <- factor(sig_g$Genus_clean, levels = tree_order_g)

phylum_order <- sig_g %>%
  arrange(match(Genus_clean, tree_order_g)) %>%
  pull(Phylum) %>%
  unique()

sig_g$Phylum <- factor(sig_g$Phylum, levels = phylum_order)
tip_tax_g$Phylum <- factor(tip_tax_g$Phylum, levels = phylum_order)

sig_g$log10baseMean <- log10(sig_g$baseMean + 1)

# Save significant genus results
sig_g_save <- sig_g
sig_g_save$Genus_clean <- as.character(sig_g_save$Genus_clean)
sig_g_save$Phylum <- as.character(sig_g_save$Phylum)

write.xlsx(sig_g_save,file.path(out_dir, "DESeq2_high_vs_no_fertilizer_bacterial_genus_Leonard_2022_metagenomics.xlsx"),
  overwrite = TRUE)

# Tree plot
p_tree_final <- ggtree(
  phylo_tree_g,
  aes(color = Phylum),
  size = 9
) %<+% tip_tax_g +
  scale_color_manual(values = phylum_colors, na.value = "grey50") +
  theme_tree() +
  theme(legend.position = "none")

# Log2 fold-change plot
p_bars_final <- ggplot(
  sig_g,
  aes(x = log2FoldChange, y = Genus_clean, fill = Phylum)
) +
  geom_bar(stat = "identity", color = "black", linewidth = 5, width = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 5) +
  scale_fill_manual(values = phylum_colors) +
  scale_x_continuous(breaks = seq(-1.5, 2, 0.5)) +
  labs(
    x = "log2FoldChange\n← No fertilizer          High fertilizer →",
    y = ""
  ) +
  theme_bw(base_size = 200) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 180, color = "black", face = "bold"),
    axis.text.x = element_text(size = 400, color = "black", face = "bold", angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 2.5),
    plot.margin = margin(10, 10, 150, 10)
  )

# BaseMean plot
p_basemean <- ggplot(
  sig_g,
  aes(x = log10baseMean, y = Genus_clean)
) +
  geom_bar(stat = "identity", fill = "#0097A7", color = "black", linewidth = 2, width = 1) +
  theme_bw(base_size = 200) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 400, color = "black", face = "bold", angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1.5),
    plot.margin = margin(10, 10, 150, 10)
  )

# Combine plots
final_genus <- p_tree_final + p_bars_final + p_basemean +
  plot_layout(widths = c(1, 3, 1.2))

plot_height <- nrow(sig_g) * 2.2


# Save figure as PDF
ggsave(file.path(
  out_dir, "DESeq2_high_vs_no_fertilizer_bacterial_genus_Leonard_2022_metagenomics.pdf"),
  final_genus,
  width = 180,
  height = plot_height,
  units = "in",
  bg = "white",
  limitsize = FALSE
)

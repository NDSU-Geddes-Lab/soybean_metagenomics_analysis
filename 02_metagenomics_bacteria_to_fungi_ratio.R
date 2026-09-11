# Bacteria-to-fungi ratio analysis across different fertilizer treatments
# in the soybean rhizosphere at Leonard in 2022 using metagenomics data

# Load required packages
library(phyloseq)
library(forcats)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(rstatix)
library(dplyr)
library(writexl)

# Output directory
out_dir <- "Bacteria_to_fungi_ratio"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Fertilizer order and colors
fert_order <- c("No fertilizer", "Medium fertilizer", "High fertilizer")

fert_box_colors <- c(
  "No fertilizer" = "#8B1C62",
  "Medium fertilizer" = "#00008B",
  "High fertilizer" = "#009E73"
)

fert_point_colors <- c(
  "No fertilizer" = "#5a0f3f",
  "Medium fertilizer" = "#00004f",
  "High fertilizer" = "#006b47"
)

# Extract bacterial read counts
ps_metag_b <- readRDS("PS.metagenomics_bacteria.rds")

bac_reads <- data.frame(
  SampleID = sample_names(ps_metag_b),
  BacteriaReads = sample_sums(ps_metag_b)
)

# Extract fungal read counts
ps_metag_f <- readRDS("PS.metagenomics_fungi.rds")

fung_reads <- data.frame(
  SampleID = sample_names(ps_metag_f),
  FungiReads = sample_sums(ps_metag_f)
)

# Merge read counts and calculate bacteria-to-fungi ratio
bf_reads <- full_join(bac_reads, fung_reads, by = "SampleID") %>%
  mutate(bf_ratio = BacteriaReads / FungiReads)

# Add metadata
meta <- as.data.frame(sample_data(ps_metag_b))
meta$SampleID <- rownames(meta)

df <- bf_reads %>%
  left_join(meta, by = "SampleID") %>%
  mutate(Fertilizer = fct_relevel(factor(Fertilizer), fert_order))

# Save bacteria-to-fungi ratio table
write_xlsx(
  df,
  file.path(out_dir, "Bacteria_to_Fungi_ratio.xlsx")
)

# Shapiro-Wilk normality tests
sw_bacteria <- do.call(rbind, lapply(fert_order, function(g) {
  
  vals <- df$BacteriaReads[df$Fertilizer == g]
  sw <- shapiro.test(vals)
  
  data.frame(
    Variable = "BacteriaReads",
    Fertilizer = g,
    n = length(vals),
    W = round(sw$statistic, 4),
    p = round(sw$p.value, 4),
    Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
  
}))

sw_fungi <- do.call(rbind, lapply(fert_order, function(g) {
  
  vals <- df$FungiReads[df$Fertilizer == g]
  sw <- shapiro.test(vals)
  
  data.frame(
    Variable = "FungiReads",
    Fertilizer = g,
    n = length(vals),
    W = round(sw$statistic, 4),
    p = round(sw$p.value, 4),
    Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
  
}))

sw_ratio <- do.call(rbind, lapply(fert_order, function(g) {
  
  vals <- df$bf_ratio[df$Fertilizer == g]
  sw <- shapiro.test(vals)
  
  data.frame(
    Variable = "bf_ratio",
    Fertilizer = g,
    n = length(vals),
    W = round(sw$statistic, 4),
    p = round(sw$p.value, 4),
    Normal = ifelse(sw$p.value > 0.05, "YES", "NO")
  )
  
}))

sw_all <- rbind(sw_bacteria, sw_fungi, sw_ratio)
print(sw_all)

# BacteriaReads - Normality not met - Kruskal-Wallis + pairwise Wilcoxon with BH correction
kw_bacteria <- df %>%
  kruskal_test(BacteriaReads ~ Fertilizer)

print(as.data.frame(kw_bacteria))

pwc_bacteria_all <- df %>%
  wilcox_test(
    BacteriaReads ~ Fertilizer,
    p.adjust.method = "BH",
    paired = FALSE
  ) %>%
  add_significance("p.adj")

print(as.data.frame(pwc_bacteria_all))

kw_bacteria_export <- as.data.frame(kw_bacteria) %>%
  mutate(
    variable = "BacteriaReads",
    test = "Kruskal-Wallis",
    group1 = NA_character_,
    group2 = NA_character_,
    p.adj = NA_real_,
    p.adj.signif = NA_character_
  ) %>%
  select(
    variable,
    test,
    group1,
    group2,
    statistic,
    p,
    p.adj,
    p.adj.signif
  )

pwc_bacteria_export <- as.data.frame(pwc_bacteria_all) %>%
  mutate(
    variable = "BacteriaReads",
    test = "Pairwise Wilcoxon (BH adjusted)"
  ) %>%
  select(
    variable,
    test,
    group1,
    group2,
    statistic,
    p,
    p.adj,
    p.adj.signif
  )

# FungiReads - Normality not met - Kruskal-Wallis + pairwise Wilcoxon with BH correction
kw_fungi <- df %>%
  kruskal_test(FungiReads ~ Fertilizer)

print(as.data.frame(kw_fungi))

pwc_fungi_all <- df %>%
  wilcox_test(
    FungiReads ~ Fertilizer,
    p.adjust.method = "BH",
    paired = FALSE
  ) %>%
  add_significance("p.adj")

print(as.data.frame(pwc_fungi_all))

kw_fungi_export <- as.data.frame(kw_fungi) %>%
  mutate(
    variable = "FungiReads",
    test = "Kruskal-Wallis",
    group1 = NA_character_,
    group2 = NA_character_,
    p.adj = NA_real_,
    p.adj.signif = NA_character_
  ) %>%
  select(
    variable,
    test,
    group1,
    group2,
    statistic,
    p,
    p.adj,
    p.adj.signif
  )

pwc_fungi_export <- as.data.frame(pwc_fungi_all) %>%
  mutate(
    variable = "FungiReads",
    test = "Pairwise Wilcoxon (BH adjusted)"
  ) %>%
  select(
    variable,
    test,
    group1,
    group2,
    statistic,
    p,
    p.adj,
    p.adj.signif
  )

# Bacteria-to-fungi ratio - Normality met - one-way ANOVA + Tukey HSD
aov_ratio <- aov(
  bf_ratio ~ Fertilizer,
  data = df
)

print(summary(aov_ratio))

tukey_ratio_all <- df %>%
  tukey_hsd(bf_ratio ~ Fertilizer) %>%
  add_significance("p.adj")

print(as.data.frame(tukey_ratio_all))

aov_ratio_export <- data.frame(
  variable = "bf_ratio",
  test = "One-way ANOVA",
  group1 = NA_character_,
  group2 = NA_character_,
  statistic = summary(aov_ratio)[[1]]$`F value`[1],
  p = summary(aov_ratio)[[1]]$`Pr(>F)`[1],
  p.adj = NA_real_,
  p.adj.signif = NA_character_
)

tukey_ratio_export <- as.data.frame(tukey_ratio_all) %>%
  mutate(
    variable = "bf_ratio",
    test = "Tukey HSD",
    statistic = NA_real_,
    p = NA_real_
  ) %>%
  select(
    variable,
    test,
    group1,
    group2,
    statistic,
    p,
    p.adj,
    p.adj.signif
  )

# Save statistical results
inferential_statistics <- bind_rows(
  kw_bacteria_export,
  pwc_bacteria_export,
  kw_fungi_export,
  pwc_fungi_export,
  aov_ratio_export,
  tukey_ratio_export
)

write_xlsx(
  list(
    Shapiro_Wilk = sw_all,
    Statistical_tests = inferential_statistics
  ),
  file.path(out_dir, "Bacteria_to_Fungi_ratio_statistics.xlsx")
)

# Position significance bracket - No fertilizer vs High fertilizer only
tukey_ratio <- tukey_ratio_all %>%
  filter(
    group1 == "No fertilizer",
    group2 == "High fertilizer",
    p.adj <= 0.05
  ) %>%
  mutate(
    p.adj.signif = ifelse(
      p.adj.signif == "****",
      "***",
      p.adj.signif
    )
  )

y_max <- max(df$bf_ratio, na.rm = TRUE)

tukey_ratio <- tukey_ratio %>%
  mutate(
    y.position = y_max + 0.25
  )

y_upper <- if (nrow(tukey_ratio) > 0) {
  max(tukey_ratio$y.position) + 0.5
} else {
  y_max + 1
}

# Figure 
p_ratio <- ggplot(
  df,
  aes(
    x = Fertilizer,
    y = bf_ratio,
    fill = Fertilizer
  )
) +
  geom_boxplot(
    alpha = 0.7,
    outlier.shape = NA,
    width = 0.6,
    color = "black",
    linewidth = 0.8
  ) +
  geom_jitter(
    aes(color = Fertilizer),
    width = 0.15,
    size = 18,
    alpha = 0.8
  ) +
  {
    if (nrow(tukey_ratio) > 0)
      stat_pvalue_manual(
        tukey_ratio,
        label = "p.adj.signif",
        y.position = "y.position",
        tip.length = 0.02,
        bracket.size = 8,
        size = 56,
        fontface = "bold",
        color = "black",
        inherit.aes = FALSE
      )
  } +
  scale_x_discrete(
    limits = fert_order,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, y_upper),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_fill_manual(
    values = fert_box_colors,
    limits = fert_order,
    drop = FALSE
  ) +
  scale_color_manual(
    values = fert_point_colors,
    limits = fert_order,
    drop = FALSE
  ) +
  labs(
    x = "",
    y = "Bacteria/Fungi Read Ratio"
  ) +
  theme_cowplot() +
  theme(
    axis.line = element_line(linewidth = 3),
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      size = 120
    ),
    axis.text.y = element_text(
      face = "bold",
      size = 120
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 120
    )
  )

# Save figure
ggsave(
  filename = file.path(out_dir, "Bacteria_to_Fungi_ratio.png"),
  plot = p_ratio,
  device = "png",
  width = 26,
  height = 33,
  units = "in",
  dpi = 600,
  bg = "white"
)

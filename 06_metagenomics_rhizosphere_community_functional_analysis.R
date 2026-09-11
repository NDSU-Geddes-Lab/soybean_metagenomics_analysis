# Differential abundance analysis using functional modules in soybean rhizosphere community 
# between high and no fertilizer treatments at Leonard in 2022


# Load required packages
library(DESeq2)
library(readxl)
library(dplyr)
library(ggplot2)
library(openxlsx)


# Output directory
out_dir <- "metagenomics_rhizosphere_community_functional_analysis"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Load metadata
meta <- read.csv(
  "metagenomics_metadata.csv",
  stringsAsFactors = FALSE
)

# Keep No fertilizer and High fertilizer treatments
meta <- meta %>%
  filter(
    Fertilizer %in% c("No fertilizer", "High fertilizer")
  )

# Match SampleID format to DRAM sample columns
meta$colname <- paste0(
  meta$SampleID,
  "_contigs"
)

# Set No fertilizer as the reference level
meta$Fertilizer <- factor(
  meta$Fertilizer,
  levels = c("No fertilizer", "High fertilizer")
)

# Load DRAM gene-level counts
sheets <- c(
  "MISC",
  "carbon utilization",
  "Transporters",
  "Energy",
  "Organic Nitrogen",
  "carbon utilization (Woodcroft)"
)

all_genes <- lapply(
  sheets,
  function(sh) {
    
    df <- read_excel(
      "metabolism_summary.xlsx",
      sheet = sh
    )
    
    df$sheet <- sh
    
    df
  }
) %>%
  bind_rows()

# Identify sample columns
sample_cols <- intersect(
  meta$colname,
  colnames(all_genes)
)

# Convert sample counts to numeric and replace missing values with zero
all_genes[sample_cols] <- lapply(
  all_genes[sample_cols],
  function(x) as.numeric(as.character(x))
)

all_genes[sample_cols][
  is.na(all_genes[sample_cols])
] <- 0

# Aggregate gene counts to functional modules
module_counts <- all_genes %>%
  group_by(
    module,
    sheet
  ) %>%
  summarise(
    across(
      all_of(sample_cols),
      sum
    ),
    .groups = "drop"
  ) %>%
  mutate(
    module_id = paste0(
      module,
      " [",
      sheet,
      "]"
    )
  )

# Keep module information
module_info <- module_counts %>%
  select(
    module_id,
    module,
    sheet
  )

# Create count matrix
count_matrix <- as.matrix(
  module_counts %>%
    select(
      all_of(sample_cols)
    )
)

rownames(count_matrix) <- module_counts$module_id

storage.mode(count_matrix) <- "integer"

stopifnot(
  !any(
    duplicated(
      rownames(count_matrix)
    )
  )
)

# Filter low-count modules
keep <- rowSums(count_matrix) >= 10

count_matrix <- count_matrix[
  keep,
]

module_info <- module_info[
  keep,
]

# Match metadata to count matrix
coldata <- meta[
  match(
    colnames(count_matrix),
    meta$colname
  ),
]

stopifnot(
  all(
    coldata$colname == colnames(count_matrix)
  )
)

# DESeq2 functional module analysis
dds_modules <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ Fertilizer
)

dds_modules <- DESeq(
  dds_modules
)

alpha <- 0.05

res_modules <- results(
  dds_modules,
  contrast = c(
    "Fertilizer",
    "High fertilizer",
    "No fertilizer"
  ),
  alpha = alpha
)

# Add functional module information
res_modules_df <- as.data.frame(
  res_modules
)

res_modules_df$module_id <- rownames(
  res_modules_df
)

res_modules_df <- res_modules_df %>%
  left_join(
    module_info,
    by = "module_id"
  ) %>%
  mutate(
    Direction = case_when(
      padj < alpha & log2FoldChange > 0 ~ "Higher in High fertilizer",
      padj < alpha & log2FoldChange < 0 ~ "Higher in No fertilizer",
      TRUE ~ "Not significant"
    )
  ) %>%
  arrange(
    padj
  )

# Keep significant functional modules
sig_modules <- res_modules_df %>%
  filter(
    !is.na(padj),
    padj < alpha
  )

# Save DESeq2 results

write.xlsx(
  list(
    All_modules = res_modules_df,
    Significant_modules = sig_modules
  ),
  file.path(
    out_dir,
    "DESeq2_high_vs_no_fertilizer_rhizosphere_functional_modules_Leonard_2022_metagenomics.xlsx"
  ),
  overwrite = TRUE
)

# Prepare significant modules for plotting
plot_df <- sig_modules %>%
  mutate(
    module_short = ifelse(
      nchar(module) > 55,
      paste0(substr(module, 1, 53), "..."),
      module
    )
  ) %>%
  arrange(
    sheet,
    log2FoldChange
  ) %>%
  mutate(
    module_short = factor(
      module_short,
      levels = unique(module_short)
    )
  )

# Functional category colors
sheet_colors <- c(
  "MISC" = "#9467bd",
  "carbon utilization" = "#2ca02c",
  "carbon utilization (Woodcroft)" = "#bcbd22",
  "Transporters" = "#1f77b4",
  "Energy" = "#e377c2",
  "Organic Nitrogen" = "#ff7f0e"
)

# Log2 fold-change plot
p_bars <- ggplot(
  plot_df,
  aes(
    x = log2FoldChange,
    y = module_short,
    fill = sheet
  )
) +
  geom_col(
    alpha = 0.9,
    width = 0.75
  ) +
  geom_vline(
    xintercept = 0,
    color = "black",
    linewidth = 0.3
  ) +
  facet_grid(
    sheet ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_fill_manual(
    values = sheet_colors,
    guide = "none"
  ) +
  labs(
    x = "← Higher in No Fertilizer    log2 Fold Change    Higher in High Fertilizer →",
    y = NULL
  ) +
  theme_minimal(
    base_size = 10
  ) +
  theme(
    axis.text.y = element_text(
      size = 7
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(
      color = "grey85"
    ),
    panel.spacing.y = grid::unit(
      0.15,
      "lines"
    ),
    strip.text.y = element_blank(),
    strip.background.y = element_blank()
  )

# Save figure
ggsave(
  filename = file.path(
    out_dir,
    "DESeq2_high_vs_no_fertilizer_functional_modules_Leonard_2022_metagenomics.png"
  ),
  plot = p_bars,
  width = 16,
  height = 11,
  units = "in",
  dpi = 600,
  bg = "white"
)

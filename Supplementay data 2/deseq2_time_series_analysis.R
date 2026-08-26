# ============================================================
# DESeq2 Analysis of Temporal Dynamics in Extracellular Enzyme Genes
#
# This workflow identifies genes exhibiting significant temporal
# changes across the incubation period and evaluates differential
# abundance between Day 0 and Day 720.
#
# The analysis consists of:
#   (i) likelihood ratio test (LRT) for overall temporal trends,
#   (ii) Wald test for pairwise comparison (Day 720 vs Day 0),
#   (iii) integration of both analyses to identify robust signals,
#   (iv) estimation of normalized abundance ratios.
#
# Time points:
#   1 = Day 0
#   2 = Day 9
#   3 = Day 70
#   4 = Day 450
#   5 = Day 720
# ============================================================

setwd("D:/ExperimentalRecords/2025_RPOC_Revision/DESeq2_Analysis")

# ---- Load required package ----

library(DESeq2)

# ============================================================
# 1. Import count matrix and sample metadata
# ============================================================

count_matrix <- read.table(
  "counts.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

sample_metadata <- read.table(
  "coldata.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

# Ensure identical sample order between matrices

count_matrix <- count_matrix[, rownames(sample_metadata)]

# Encode time as a continuous variable

sample_metadata$time <- as.numeric(sample_metadata$time)

cat(
  "Time levels detected:",
  unique(sample_metadata$time),
  "\n"
)

# Expected coding:
# 1 = Day 0
# 2 = Day 9
# 3 = Day 70
# 4 = Day 450
# 5 = Day 720

# ============================================================
# 2. Likelihood ratio test (LRT)
# Identify genes with significant temporal trends
# ============================================================

dds_lrt <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_metadata,
  design = ~ time
)

# Remove low-abundance genes

dds_lrt <- dds_lrt[
  rowSums(counts(dds_lrt)) > 5,
]

dds_lrt <- DESeq(
  dds_lrt,
  test = "LRT",
  reduced = ~ 1
)

# Extract regression coefficients for time

lrt_results <- results(
  dds_lrt,
  name = "time"
)

lrt_df <- as.data.frame(lrt_results)

lrt_df$gene <- rownames(lrt_df)

lrt_df$trend <- ifelse(
  lrt_df$log2FoldChange > 0,
  "Upregulated",
  ifelse(
    lrt_df$log2FoldChange < 0,
    "Downregulated",
    "No_change"
  )
)

lrt_df <- lrt_df[!is.na(lrt_df$padj), ]
lrt_df <- lrt_df[order(lrt_df$padj), ]

significant_lrt <- subset(
  lrt_df,
  padj < 0.05
)

write.csv(
  lrt_df,
  "DESeq2_temporal_trends_all_genes.csv",
  row.names = FALSE
)

write.csv(
  significant_lrt,
  "DESeq2_temporal_trends_significant_genes.csv",
  row.names = FALSE
)

cat(
  "LRT completed: significant genes =",
  nrow(significant_lrt),
  "\n"
)

# ============================================================
# 3. Pairwise comparison: Day 720 versus Day 0
# ============================================================

selected_samples <- rownames(sample_metadata)[
  sample_metadata$time %in% c(1, 5)
]

count_matrix_pair <- count_matrix[, selected_samples]

metadata_pair <- sample_metadata[
  selected_samples,
  ,
  drop = FALSE
]

cat(
  "Number of retained samples:",
  length(selected_samples),
  "\n"
)

print(table(metadata_pair$time))

# Convert time to a categorical factor

metadata_pair$time_factor <- factor(
  metadata_pair$time
)

metadata_pair$time_factor <- relevel(
  metadata_pair$time_factor,
  ref = "1"
)

cat(
  "Reference level:",
  levels(metadata_pair$time_factor)[1],
  "\n"
)

dds_pair <- DESeqDataSetFromMatrix(
  countData = count_matrix_pair,
  colData = metadata_pair,
  design = ~ time_factor
)

dds_pair <- dds_pair[
  rowSums(counts(dds_pair)) > 5,
]

dds_pair <- DESeq(dds_pair)

cat("Available coefficients:\n")
print(resultsNames(dds_pair))

pairwise_results <- results(
  dds_pair,
  contrast = c(
    "time_factor",
    "5",
    "1"
  ),
  alpha = 0.05
)

pairwise_df <- as.data.frame(pairwise_results)

pairwise_df$gene <- rownames(pairwise_df)

pairwise_df$trend <- ifelse(
  pairwise_df$log2FoldChange > 0,
  "Upregulated",
  ifelse(
    pairwise_df$log2FoldChange < 0,
    "Downregulated",
    "No_change"
  )
)

pairwise_df <- pairwise_df[
  !is.na(pairwise_df$padj),
]

pairwise_df <- pairwise_df[
  order(pairwise_df$padj),
]

significant_pairwise <- subset(
  pairwise_df,
  padj < 0.05
)

write.csv(
  pairwise_df,
  "DESeq2_Day720_vs_Day0_all_genes.csv",
  row.names = FALSE
)

write.csv(
  significant_pairwise,
  "DESeq2_Day720_vs_Day0_significant_genes.csv",
  row.names = FALSE
)

cat(
  "Pairwise comparison completed: significant genes =",
  nrow(significant_pairwise),
  "\n"
)

# ============================================================
# 4. Integrate LRT and pairwise results
# Retain genes with consistent directional changes
# ============================================================

integrated_results <- merge(
  significant_lrt[
    ,
    c(
      "gene",
      "log2FoldChange",
      "padj",
      "trend"
    )
  ],
  significant_pairwise[
    ,
    c(
      "gene",
      "log2FoldChange",
      "padj",
      "trend"
    )
  ],
  by = "gene",
  suffixes = c(
    "_lrt",
    "_pairwise"
  )
)

consistent_genes <- subset(
  integrated_results,
  trend_lrt == trend_pairwise
)

consistent_genes <- consistent_genes[
  order(consistent_genes$padj_pairwise),
]

write.csv(
  consistent_genes,
  "DESeq2_consistent_temporal_genes.csv",
  row.names = FALSE
)

cat(
  "Consistently significant genes:",
  nrow(consistent_genes),
  "\n"
)

cat(
  "Upregulated genes:",
  sum(consistent_genes$trend_lrt == "Upregulated"),
  "\n"
)

cat(
  "Downregulated genes:",
  sum(consistent_genes$trend_lrt == "Downregulated"),
  "\n"
)

# ============================================================
# 5. Calculate normalized abundance ratios
#
# Note:
# DESeq2 normalized counts are not TPM values.
# Fold changes below are calculated using size factor-
# normalized counts and are intended for descriptive purposes.
# ============================================================

normalized_counts <- counts(
  dds_lrt,
  normalized = TRUE
)

day0_samples <- rownames(sample_metadata)[
  sample_metadata$time == 1
]

day720_samples <- rownames(sample_metadata)[
  sample_metadata$time == 5
]

mean_day0 <- rowMeans(
  normalized_counts[
    ,
    day0_samples,
    drop = FALSE
  ]
)

mean_day720 <- rowMeans(
  normalized_counts[
    ,
    day720_samples,
    drop = FALSE
  ]
)

abundance_summary <- data.frame(
  gene = rownames(normalized_counts),

  mean_day0 = round(mean_day0, 2),

  mean_day720 = round(mean_day720, 2),

  normalized_fold_change = round(
    (mean_day720 + 1) /
      (mean_day0 + 1),
    2
  )
)

final_results <- merge(
  consistent_genes,
  abundance_summary,
  by = "gene"
)

final_results <- final_results[
  order(final_results$padj_pairwise),
]

write.csv(
  final_results,
  "DESeq2_final_results.csv",
  row.names = FALSE
)

# ============================================================
# 6. Summary of output files
# ============================================================

cat("\nAnalysis completed successfully.\n\n")

cat("Generated files:\n")

cat(
  "  DESeq2_temporal_trends_all_genes.csv\n"
)

cat(
  "  DESeq2_temporal_trends_significant_genes.csv\n"
)

cat(
  "  DESeq2_Day720_vs_Day0_all_genes.csv\n"
)

cat(
  "  DESeq2_Day720_vs_Day0_significant_genes.csv\n"
)

cat(
  "  DESeq2_consistent_temporal_genes.csv\n"
)

cat(
  "  DESeq2_final_results.csv\n"
)
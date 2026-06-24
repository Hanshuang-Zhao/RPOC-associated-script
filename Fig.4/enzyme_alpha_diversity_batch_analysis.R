# ============================================================
# Alpha Diversity Analysis of CAZyme-, Peptidase-, and
# Microbial Community Profiles
#
# This workflow calculates alpha diversity metrics for
# multiple feature tables, including:
#
#   1. Secreted CAZymes
#   2. Total CAZymes
#   3. Secreted Peptidases
#   4. Total Peptidases
#   5. Secreted-CAZyme-associated microbial community
#   6. Total-CAZyme-associated microbial community
#   7. Secreted-peptidase-associated microbial community
#   8. Total-peptidase-associated microbial community
#
# Metrics:
#   - Observed richness
#   - Shannon index
#   - Pielou's evenness
#   - Gini–Simpson index
#   - Inverse Simpson index
#   - Chao1 richness estimator
#   - ACE richness estimator
#   - Good's coverage
#
# ============================================================

setwd("D:/ExperimentalRecords/2025_RPOC_Revision/AlphaDiversity")

library(vegan)

# ============================================================
# Input files
# ============================================================

input_files <- c(
  "secreted_CAZymes.txt",
  "total_CAZymes.txt",
  "secreted_peptidases.txt",
  "total_peptidases.txt",
  "secreted_CAZyme_microbiome.txt",
  "total_CAZyme_microbiome.txt",
  "secreted_peptidase_microbiome.txt",
  "total_peptidase_microbiome.txt"
)

# ============================================================
# Alpha diversity function
# ============================================================

calculate_alpha_diversity <- function(file_name){

  cat("\nProcessing:", file_name, "\n")

  feature_table <- read.delim(
    file_name,
    row.names = 1,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # samples × features

  feature_table <- t(feature_table)
  feature_table <- as.matrix(feature_table)

  # ----------------------------------------------------------
  # Richness
  # ----------------------------------------------------------

  richness <- rowSums(feature_table > 0)

  # ----------------------------------------------------------
  # Shannon
  # ----------------------------------------------------------

  shannon_index <- diversity(
    feature_table,
    index = "shannon"
  )

  effective_shannon <- exp(shannon_index)

  pielou_evenness <- shannon_index /
    log(richness)

  # ----------------------------------------------------------
  # Simpson
  # ----------------------------------------------------------

  gini_simpson <- diversity(
    feature_table,
    index = "simpson"
  )

  inverse_simpson <- diversity(
    feature_table,
    index = "invsimpson"
  )

  simpson_equitability <- inverse_simpson /
    richness

  # ----------------------------------------------------------
  # Chao1 and ACE
  # ----------------------------------------------------------

  feature_table_int <- round(feature_table)

  richness_est <- estimateR(
    feature_table_int
  )

  chao1 <- richness_est[2, ]

  ace <- richness_est[4, ]

  # ----------------------------------------------------------
  # Good's coverage
  # ----------------------------------------------------------

  goods_coverage <-
    1 -
    rowSums(feature_table_int == 1) /
    rowSums(feature_table_int)

  # ----------------------------------------------------------
  # Output table
  # ----------------------------------------------------------

  alpha_results <- data.frame(
    sample_id = rownames(feature_table),

    richness = richness,

    shannon_index = shannon_index,

    effective_shannon = effective_shannon,

    pielou_evenness = pielou_evenness,

    gini_simpson = gini_simpson,

    inverse_simpson = inverse_simpson,

    simpson_equitability =
      simpson_equitability,

    chao1 = chao1,

    ace = ace,

    goods_coverage =
      goods_coverage
  )

  output_name <- paste0(
    tools::file_path_sans_ext(file_name),
    "_alpha_diversity.csv"
  )

  write.csv(
    alpha_results,
    output_name,
    row.names = FALSE
  )

  cat(
    "Output saved:",
    output_name,
    "\n"
  )

  return(alpha_results)
}

# ============================================================
# Run all analyses
# ============================================================

alpha_diversity_results <- lapply(
  input_files,
  calculate_alpha_diversity
)

names(alpha_diversity_results) <-
  tools::file_path_sans_ext(input_files)

cat(
  "\nAll alpha diversity analyses completed.\n"
)
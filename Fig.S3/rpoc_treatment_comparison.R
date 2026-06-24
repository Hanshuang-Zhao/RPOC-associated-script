# ============================================================
# Statistical Comparison of RPOC Concentrations Between
# Control and Inoculated Treatments
#
# This script performs group-specific comparisons of RPOC
# concentrations using:
#
#   (i) Welch's two-sample t-test
#   (ii) Mann–Whitney U test (Wilcoxon rank-sum test)
#
# Welch's t-test was selected because it does not assume
# equal variances between treatments.
#
# The non-parametric Mann–Whitney U test was additionally
# conducted to evaluate the robustness of the results when
# normality assumptions may not be satisfied.
# ============================================================

setwd("D:/ExperimentalRecords/2025_RPOC_Revision/DataProcessing")

# ---- Load required packages ----

library(dplyr)

# ============================================================
# 1. Import data and inspect structure
# ============================================================

rpoc_data <- read.table(
  "RPOC_control_data.txt",
  header = TRUE,
  sep = "\t"
)

str(rpoc_data)

table(
  rpoc_data$Group,
  rpoc_data$Treatment
)

# Ensure a consistent treatment order

rpoc_data$Treatment <- factor(
  rpoc_data$Treatment,
  levels = c("Control", "Inoculated")
)

# ============================================================
# 2. Define group-specific statistical tests
# ============================================================

run_group_comparison <- function(data, group_name) {

  subset_data <- data %>%
    filter(Group == group_name)

  control_values <- subset_data$RPOC[
    subset_data$Treatment == "Control"
  ]

  inoculated_values <- subset_data$RPOC[
    subset_data$Treatment == "Inoculated"
  ]

  # Welch's t-test

  welch_test <- t.test(
    control_values,
    inoculated_values,
    var.equal = FALSE
  )

  # Mann–Whitney U test

  wilcox_test <- wilcox.test(
    control_values,
    inoculated_values,
    exact = FALSE
  )

  data.frame(
    group = group_name,

    n_control = length(control_values),
    n_inoculated = length(inoculated_values),

    mean_control = mean(control_values, na.rm = TRUE),
    mean_inoculated = mean(inoculated_values, na.rm = TRUE),

    mean_difference =
      mean(inoculated_values, na.rm = TRUE) -
      mean(control_values, na.rm = TRUE),

    t_statistic = unname(welch_test$statistic),
    degrees_freedom = unname(welch_test$parameter),
    p_value_ttest = welch_test$p.value,

    W_statistic = unname(wilcox_test$statistic),
    p_value_wilcox = wilcox_test$p.value
  )
}

# ============================================================
# 3. Apply statistical tests to each group
# ============================================================

group_names <- unique(rpoc_data$Group)

comparison_results <- bind_rows(
  lapply(group_names, function(group_name) {
    run_group_comparison(
      rpoc_data,
      group_name
    )
  })
)

# ============================================================
# 4. Determine statistical significance
# ============================================================

comparison_results <- comparison_results %>%
  mutate(
    significant_ttest = ifelse(
      p_value_ttest < 0.05,
      "*",
      "ns"
    ),

    significant_wilcox = ifelse(
      p_value_wilcox < 0.05,
      "*",
      "ns"
    )
  )

# ============================================================
# 5. Display results
# ============================================================

cat("\n====================================================\n")
cat("Group-specific comparisons of RPOC concentrations\n")
cat("====================================================\n")

print(comparison_results)

# ============================================================
# 6. Export results
# ============================================================

write.table(
  comparison_results,
  file = "RPOC_group_comparison_results.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nAnalysis completed successfully.\n")
cat("Results saved as: RPOC_group_comparison_results.txt\n")
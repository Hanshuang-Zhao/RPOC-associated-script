# ============================================================
# Biphasic POC Degradation Analysis Across Experimental Treatments
#
# This script quantifies particulate organic carbon (POC) degradation
# dynamics using a two-phase exponential decay framework.
#
# Analyses include:
#   (i) estimation of phase-specific decay constants (k),
#   (ii) within-treatment comparisons between degradation phases,
#   (iii) between-treatment comparisons within each phase, and
#   (iv) pairwise ANCOVA tests of slope heterogeneity.
#
# Statistical procedures follow the approach described in:
# Craft & Bochdansky (2025).
# ============================================================

setwd("D:/ExperimentalRecords/2025_RPOC_Revision/DataProcessing")

# ---- Load required packages ----

library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# 1. Import and preprocess POC datasets
# ============================================================

read_poc_data <- function(file_path, treatment_label) {

  poc_data <- read.table(
    file_path,
    header = TRUE,
    sep = "\t"
  )

  poc_long <- poc_data %>%
    pivot_longer(
      cols = -Time,
      names_to = "replicate",
      values_to = "POC"
    ) %>%
    group_by(replicate) %>%
    mutate(
      normalized_poc = POC / POC[Time == 0],
      percent_remaining = normalized_poc * 100,
      ln_percent_remaining = log(percent_remaining)
    ) %>%
    ungroup() %>%
    mutate(
      treatment = treatment_label,
      phase = case_when(
        Time <= 20 ~ "Phase1",
        Time >= 50 ~ "Phase2",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(phase))

  return(poc_long)
}

long_term_coculture <- read_poc_data(
  "POC_longterm.txt",
  "720d_coculture"
)

short_term_coculture <- read_poc_data(
  "POC_shortterm.txt",
  "30d_coculture"
)

axenic_culture <- read_poc_data(
  "POC_axenic.txt",
  "Axenic"
)

all_data <- bind_rows(
  long_term_coculture,
  short_term_coculture,
  axenic_culture
)

# Define factor order for plotting and statistical comparisons

all_data$treatment <- factor(
  all_data$treatment,
  levels = c(
    "30d_coculture",
    "720d_coculture",
    "Axenic"
  )
)

all_data$phase <- factor(
  all_data$phase,
  levels = c("Phase1", "Phase2")
)

# ============================================================
# 2. Estimate phase-specific decay constants
# ============================================================

estimate_decay_parameters <- function(data, treatment_name, phase_name) {

  subset_data <- data %>%
    filter(
      treatment == treatment_name,
      phase == phase_name
    )

  if (nrow(subset_data) < 3) return(NULL)

  model <- lm(
    ln_percent_remaining ~ Time,
    data = subset_data
  )

  model_summary <- summary(model)

  data.frame(
    treatment = treatment_name,
    phase = phase_name,
    k = coef(model)[2],
    k_abs = abs(coef(model)[2]),
    r_squared = model_summary$r.squared,
    p_value = model_summary$coefficients[2, 4]
  )
}

treatments <- c(
  "720d_coculture",
  "30d_coculture",
  "Axenic"
)

phases <- c(
  "Phase1",
  "Phase2"
)

decay_summary <- do.call(
  rbind,
  lapply(treatments, function(treatment_name) {
    do.call(
      rbind,
      lapply(phases, function(phase_name) {
        estimate_decay_parameters(
          all_data,
          treatment_name,
          phase_name
        )
      })
    )
  })
)

rownames(decay_summary) <- NULL

cat("\n====================================================\n")
cat("Phase-specific POC decay constants\n")
cat("====================================================\n")

print(decay_summary)

# ============================================================
# 3. Within-treatment comparison:
#    Phase 1 versus Phase 2
# ============================================================

cat("\n====================================================\n")
cat("ANCOVA: comparison of degradation phases within treatments\n")
cat("====================================================\n")

for (treatment_name in treatments) {

  subset_data <- all_data %>%
    filter(treatment == treatment_name)

  ancova_model <- lm(
    ln_percent_remaining ~ Time * phase,
    data = subset_data
  )

  cat("\nTreatment:", treatment_name, "\n")

  print(anova(ancova_model))
}

# ============================================================
# 4. Between-treatment comparison within each phase
# ============================================================

cat("\n====================================================\n")
cat("Phase 1: ANCOVA homogeneity-of-slopes test\n")
cat("====================================================\n")

phase1_data <- all_data %>%
  filter(phase == "Phase1")

phase1_model <- lm(
  ln_percent_remaining ~ Time * treatment,
  data = phase1_data
)

print(anova(phase1_model))

cat("\n====================================================\n")
cat("Phase 2: ANCOVA homogeneity-of-slopes test\n")
cat("====================================================\n")

phase2_data <- all_data %>%
  filter(phase == "Phase2")

phase2_model <- lm(
  ln_percent_remaining ~ Time * treatment,
  data = phase2_data
)

print(anova(phase2_model))

# ============================================================
# 5. Pairwise ANCOVA comparisons
#
# Multiple-testing correction was not applied to maintain
# consistency with Craft & Bochdansky (2025).
# ============================================================

treatment_pairs <- list(
  c("30d_coculture", "720d_coculture"),
  c("30d_coculture", "Axenic"),
  c("720d_coculture", "Axenic")
)

run_pairwise_ancova <- function(data, phase_name, treatment_pair) {

  subset_data <- data %>%
    filter(
      phase == phase_name,
      treatment %in% treatment_pair
    ) %>%
    mutate(
      treatment = factor(
        treatment,
        levels = treatment_pair
      )
    )

  model <- lm(
    ln_percent_remaining ~ Time * treatment,
    data = subset_data
  )

  anova_results <- anova(model)

  data.frame(
    phase = phase_name,
    comparison = paste(
      treatment_pair[1],
      "vs",
      treatment_pair[2]
    ),
    F_value = round(
      anova_results["Time:treatment", "F value"],
      3
    ),
    df1 = anova_results["Time:treatment", "Df"],
    df2 = anova_results["Residuals", "Df"],
    p_value = signif(
      anova_results["Time:treatment", "Pr(>F)"],
      4
    ),
    significant = ifelse(
      anova_results["Time:treatment", "Pr(>F)"] < 0.05,
      "*",
      "ns"
    )
  )
}

pairwise_phase1 <- do.call(
  rbind,
  lapply(treatment_pairs, function(pair) {
    run_pairwise_ancova(
      all_data,
      "Phase1",
      pair
    )
  })
)

pairwise_phase2 <- do.call(
  rbind,
  lapply(treatment_pairs, function(pair) {
    run_pairwise_ancova(
      all_data,
      "Phase2",
      pair
    )
  })
)

cat("\nPhase 1 pairwise comparisons:\n")
print(pairwise_phase1, row.names = FALSE)

cat("\nPhase 2 pairwise comparisons:\n")
print(pairwise_phase2, row.names = FALSE)

pairwise_summary <- rbind(
  pairwise_phase1,
  pairwise_phase2
)

cat("\n====================================================\n")
cat("Summary table (Supporting Information Table S2)\n")
cat("====================================================\n")

print(pairwise_summary, row.names = FALSE)

# ============================================================
# 6. Generate visualization
# ============================================================

generate_predictions <- function(
    data,
    treatment_name,
    phase_name,
    time_range) {

  subset_data <- data %>%
    filter(
      treatment == treatment_name,
      phase == phase_name
    )

  model <- lm(
    ln_percent_remaining ~ Time,
    data = subset_data
  )

  prediction_data <- data.frame(
    Time = seq(
      time_range[1],
      time_range[2],
      length.out = 100
    )
  )

  prediction_data$ln_percent_remaining <- predict(
    model,
    newdata = prediction_data
  )

  prediction_data$treatment <- treatment_name
  prediction_data$phase <- phase_name

  return(prediction_data)
}

prediction_lines <- bind_rows(

  generate_predictions(
    all_data,
    "720d_coculture",
    "Phase1",
    c(0, 20)
  ),

  generate_predictions(
    all_data,
    "720d_coculture",
    "Phase2",
    c(50, 500)
  ),

  generate_predictions(
    all_data,
    "30d_coculture",
    "Phase1",
    c(0, 20)
  ),

  generate_predictions(
    all_data,
    "30d_coculture",
    "Phase2",
    c(50, 500)
  ),

  generate_predictions(
    all_data,
    "Axenic",
    "Phase1",
    c(0, 20)
  ),

  generate_predictions(
    all_data,
    "Axenic",
    "Phase2",
    c(50, 500)
  )
)

treatment_colors <- c(
  "720d_coculture" = "#e31a1c",
  "30d_coculture" = "#1f78b4",
  "Axenic" = "#33a02c"
)

treatment_labels <- c(
  "720d_coculture" = "720-day coculture",
  "30d_coculture" = "30-day coculture",
  "Axenic" = "Axenic Synechococcus"
)

ggplot() +

  geom_point(
    data = all_data,
    aes(
      x = Time,
      y = ln_percent_remaining,
      color = treatment,
      shape = phase
    ),
    size = 3,
    alpha = 0.8
  ) +

  geom_line(
    data = prediction_lines,
    aes(
      x = Time,
      y = ln_percent_remaining,
      color = treatment,
      linetype = phase
    ),
    linewidth = 1.2
  ) +

  scale_color_manual(
    values = treatment_colors,
    labels = treatment_labels,
    name = "Treatment"
  ) +

  scale_shape_manual(
    values = c(
      "Phase1" = 16,
      "Phase2" = 17
    ),
    labels = c(
      "Phase1" = "Phase 1 (0–20 d)",
      "Phase2" = "Phase 2 (50–500 d)"
    ),
    name = "Degradation phase"
  ) +

  scale_linetype_manual(
    values = c(
      "Phase1" = "solid",
      "Phase2" = "dashed"
    ),
    labels = c(
      "Phase1" = "Phase 1",
      "Phase2" = "Phase 2"
    ),
    name = "Degradation phase"
  ) +

  labs(
    x = "Time (days)",
    y = expression(ln~"(% POC remaining)")
  ) +

  theme_bw() +

  theme(
    text = element_text(size = 14),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(
  filename = "POC_biphasic_decay_combined.pdf",
  width = 8,
  height = 5
)

cat("\nAnalysis completed successfully.\n")
cat("Figure saved as: POC_biphasic_decay_combined.pdf\n")
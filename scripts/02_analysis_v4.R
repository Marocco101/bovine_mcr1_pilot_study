# ============================================
# Co-occurrence Analysis v4
# Binomial tests using RESAPATH baseline
# Author: Makiko FUJITA-SUZANNE
# Date: 10 April 2026
# ============================================
# Baseline: RESAPATH (Bovine Digestive E. coli)
# Period: 2011-2018 (reference PRJNA1166088)
# Source: https://resapath.anses.fr/
# ============================================

combined <- read.csv("~/amr_pilot/results/combined_resfinder_results.csv",
                     stringsAsFactors = FALSE)

targets <- c("mcr-1.1", "sul1", "sul2", "dfrA1", "dfrA12", "dfrA17")

all_samples <- unique(combined$Sample)

# ============================================
# Presence/Absence Matrix
# ============================================

pa_matrix <- data.frame(Sample = all_samples, stringsAsFactors = FALSE)

for (i in seq_along(targets)) {
  gene <- targets[i]
  col_name <- make.names(gene)
  samples_with_gene <- combined$Sample[combined$Resistance.gene == gene]
  pa_matrix[[col_name]] <- as.integer(all_samples %in% samples_with_gene)
}

# --- mcr-1 all variants ---
mcr1_samples <- combined$Sample[grep("^mcr-1", combined$Resistance.gene)]
pa_matrix$mcr1_any <- as.integer(all_samples %in% mcr1_samples)

# --- Aggregated columns ---
pa_matrix$sul_any <- as.integer(pa_matrix$sul1 == 1 | pa_matrix$sul2 == 1)

pa_matrix$dfrA_any <- as.integer(
  pa_matrix$dfrA1 == 1 | pa_matrix$dfrA12 == 1 | pa_matrix$dfrA17 == 1
)

pa_matrix$trap <- as.integer(pa_matrix$sul_any == 1 & pa_matrix$dfrA_any == 1)

# ============================================
# Dataset Summary
# ============================================

cat("=== Dataset Summary ===\n")
cat("PRJNA1166088: mcr-1 positive E. coli from French cattle\n")
cat("Source: Haenni et al., 2025\n")
cat("Collection period: 2011-2019\n")
cat("N =", nrow(pa_matrix), "isolates\n")

# ============================================
# Detection Rates
# ============================================

cat("\n=== Detection Rates ===\n\n")

report_cols <- c("mcr.1.1", "sul1", "sul2", "sul_any",
                 "dfrA1", "dfrA12", "dfrA17", "dfrA_any", "trap")
report_labels <- c("mcr-1.1", "sul1", "sul2", "sul_any",
                   "dfrA1", "dfrA12", "dfrA17", "dfrA_any", "trap")

detection <- sapply(report_cols, function(col) {
  n <- sum(pa_matrix[[col]])
  pct <- round(n / nrow(pa_matrix) * 100, 1)
  c(count = n, total = nrow(pa_matrix), percent = pct)
})
colnames(detection) <- report_labels

print(detection)

# ============================================
# RESAPATH Baseline Values (2011-2019)
# Source: RESAPATH Dashboard
# Population: Bovine Digestive E. coli (Jeunes)
# URL: https://resapath.anses.fr/
# ============================================

cat("\n=== RESAPATH Baseline (Bovine E. coli, 2011-2019) ===\n\n")

# --- Sulfamides resistance rate by year ---
sul_by_year <- c(
  "2011" = 86, "2012" = 84, "2013" = 81,
  "2014" = 80, "2015" = 81, "2016" = 78,
  "2017" = 75, "2018" = 78
)

# --- Trimethoprime resistance rate by year ---
tmp_by_year <- c(
  "2011" = 22, "2012" = 28, "2013" = 14,
  "2014" = 39, "2015" = 32, "2016" = 34,
  "2017" = 35, "2018" = 37
)

# --- Trimethoprime-Sulfamides combined resistance rate by year ---
tmpsul_by_year <- c(
  "2011" = 37, "2012" = 38, "2013" = 39,
  "2014" = 38, "2015" = 37, "2016" = 37,
  "2017" = 38, "2018" = 40
)

# --- Calculate means ---
baseline_sul  <- round(mean(sul_by_year), 1)
baseline_tmp  <- round(mean(tmp_by_year), 1)
baseline_tmpsul <- round(mean(tmpsul_by_year), 1)

cat("Sulfamides:              ", baseline_sul, "%\n")
cat("Trimethoprime:           ", baseline_tmp, "%\n")
cat("Trimethoprime-Sulfamides:", baseline_tmpsul, "%\n")

# ============================================
# Fisher's Exact Test: sul vs dfrA
# ============================================

cat("\n=== Fisher's Exact Test: sul vs dfrA ===\n")
cat("Q: Do sul and dfrA co-occur more than expected?\n\n")

tbl_sul_dfr <- table(
  sul  = factor(pa_matrix$sul_any, levels = c(0, 1)),
  dfrA = factor(pa_matrix$dfrA_any, levels = c(0, 1))
)
print(tbl_sul_dfr)

fisher_sul_dfr <- fisher.test(tbl_sul_dfr)
cat("\np-value:", fisher_sul_dfr$p.value, "\n")
cat("Odds ratio:", fisher_sul_dfr$estimate, "\n")

if (fisher_sul_dfr$p.value < 0.05) {
  cat("Result: SIGNIFICANT (p < 0.05)\n")
} else {
  cat("Result: Not significant\n")
  cat("Note: sul detection rate is 97.6%, leaving insufficient\n")
  cat("variation for contingency table analysis.\n")
}

# ============================================
# Binomial Tests (RESAPATH baseline)
# ============================================

cat("\n=== Binomial Tests vs RESAPATH Baseline ===\n")
cat("Baseline: Bovine Digestive E. coli, 2011-2019 mean\n\n")

# --- Function to run and report binomial test ---
run_binom <- function(observed_col, expected_pct, label) {
  n_total <- nrow(pa_matrix)
  n_pos   <- sum(pa_matrix[[observed_col]])
  obs_pct <- round(n_pos / n_total * 100, 1)
  exp_prop <- expected_pct / 100
  
  result <- binom.test(
    x = n_pos,
    n = n_total,
    p = exp_prop,
    alternative = "greater"
  )
  
  cat("---", label, "---\n")
  cat("  Observed:", n_pos, "/", n_total, "(", obs_pct, "%)\n")
  cat("  Expected (RESAPATH):", expected_pct, "%\n")
  cat("  p-value:", result$p.value, "\n")
  
  if (result$p.value < 0.001) {
    cat("  Result: *** HIGHLY SIGNIFICANT (p < 0.001)\n\n")
  } else if (result$p.value < 0.05) {
    cat("  Result: * SIGNIFICANT (p < 0.05)\n\n")
  } else {
    cat("  Result: Not significant\n\n")
  }
  
  return(list(
    label = label,
    observed_n = n_pos,
    observed_pct = obs_pct,
    expected_pct = expected_pct,
    p.value = result$p.value,
    significant = result$p.value < 0.05
  ))
}

# Test 1: sul_any vs Sulfamides baseline
r1 <- run_binom("sul_any", baseline_sul, "sul (any) vs Sulfamides")

# Test 2: dfrA_any vs Trimethoprime baseline
r2 <- run_binom("dfrA_any", baseline_tmp, "dfrA (any) vs Trimethoprime")

# Test 3: trap vs TMP-Sulfamides baseline
r3 <- run_binom("trap", baseline_tmpsul, "Genetic Trap vs TMP-Sulfamides")

# ============================================
# Save Results
# ============================================

# --- Detection rates ---
det_df <- data.frame(
  Gene = report_labels,
  Count = as.integer(detection["count", ]),
  Total = as.integer(detection["total", ]),
  Percent = detection["percent", ],
  stringsAsFactors = FALSE
)

write.csv(det_df, "~/amr_pilot/results/detection_rates.csv",
          row.names = FALSE)

# --- Binomial test results ---
binom_summary <- data.frame(
  Test = c(r1$label, r2$label, r3$label),
  Observed_n = c(r1$observed_n, r2$observed_n, r3$observed_n),
  Observed_pct = c(r1$observed_pct, r2$observed_pct, r3$observed_pct),
  RESAPATH_baseline_pct = c(r1$expected_pct, r2$expected_pct, r3$expected_pct),
  p.value = c(r1$p.value, r2$p.value, r3$p.value),
  Significant = c(r1$significant, r2$significant, r3$significant),
  stringsAsFactors = FALSE
)

write.csv(binom_summary, "~/amr_pilot/results/binomial_tests_v4.csv",
          row.names = FALSE)

# --- Presence/absence matrix ---
write.csv(pa_matrix, "~/amr_pilot/results/presence_absence_matrix.csv",
          row.names = FALSE)

cat("=== Analysis v4 Complete ===\n")
cat("Files saved:\n")
cat("  - detection_rates.csv\n")
cat("  - binomial_tests_v4.csv\n")
cat("  - presence_absence_matrix.csv\n")
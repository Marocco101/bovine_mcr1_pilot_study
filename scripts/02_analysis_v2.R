# ============================================
# Co-occurrence Analysis v2 (Revised)
# Author: Makiko FUJITA-SUZANNE
# Date: 9 April 2026
# Note: All 127 isolates are mcr-1 positive
# ============================================

pa_matrix <- read.csv("~/amr_pilot/results/presence_absence_matrix.csv",
                      stringsAsFactors = FALSE)

cat("=== Dataset Summary ===\n")
cat("This dataset (PRJNA1166088) consists entirely of mcr-1\n")
cat("positive E. coli from French cattle (Haenni et al., 2025)\n")
cat("N =", nrow(pa_matrix), "isolates\n\n")

# ============================================
# 1. Detection Rates (Descriptive)
# ============================================

cat("=== Detection Rates among mcr-1 positive isolates ===\n\n")

genes <- c("mcr.1.1", "sul1", "sul2", "sul_any", 
           "dfrA1", "dfrA12", "dfrA17", "dfrA_any", "trap")

for (g in genes) {
  n <- sum(pa_matrix[[g]])
  pct <- round(n / nrow(pa_matrix) * 100, 1)
  cat(sprintf("  %-10s: %3d / %d  (%5.1f%%)\n", g, n, nrow(pa_matrix), pct))
}

# ============================================
# 2. Fisher's Exact Test: sul vs dfrA
# ============================================

cat("\n=== Fisher's Exact Test: sul vs dfrA ===\n")
cat("Q: Does sul+isolate tend to coexist with dfrA? \n\n")

tbl_sul_dfr <- table(
  sul  = factor(pa_matrix$sul_any, levels = c(0, 1)),
  dfrA = factor(pa_matrix$dfrA_any, levels = c(0, 1))
)
print(tbl_sul_dfr)

fisher_sul_dfr <- fisher.test(tbl_sul_dfr)
cat("\np-value:", fisher_sul_dfr$p.value, "\n")
cat("Odds ratio:", fisher_sul_dfr$estimate, "\n")

if (fisher_sul_dfr$p.value < 0.05) {
  cat("Result: SIGNIFICANT\n")
  cat("Interpretation: sul and dfrA co-occur significantly,\n")
  cat("supporting physical linkage (Genetic Trap architecture)\n")
} else {
  cat("Result: Not significant\n")
}

# ============================================
# 3. Binomial Test: 
# Compare detection rate with literature 
# ============================================

cat("\n=== Binomial Tests ===\n")
cat("Q: What is the detection rate of sul/dfrA in mcr-1 positive isolates? \n")

# Reference Detection Rate from Literature
# sul: ~40-60% in E. coli from cattle
# dfrA: ~30-50% in E. coli from cattle
expected_sul  <- 0.50
expected_dfrA <- 0.40

# Binominal Tests of sul
binom_sul <- binom.test(
  x = sum(pa_matrix$sul_any),
  n = nrow(pa_matrix),
  p = expected_sul,
  alternative = "greater"
)
cat("--- sul (expected:", expected_sul * 100, "%) ---\n")
cat("Observed:", round(mean(pa_matrix$sul_any) * 100, 1), "%\n")
cat("p-value:", binom_sul$p.value, "\n")
if (binom_sul$p.value < 0.05) {
  cat("Result: SIGNIFICANTLY HIGHER than expected\n\n")
} else {
  cat("Result: Not significantly higher\n\n")
}

# Binominal Tests of dfrA
binom_dfr <- binom.test(
  x = sum(pa_matrix$dfrA_any),
  n = nrow(pa_matrix),
  p = expected_dfrA,
  alternative = "greater"
)
cat("--- dfrA (expected:", expected_dfrA * 100, "%) ---\n")
cat("Observed:", round(mean(pa_matrix$dfrA_any) * 100, 1), "%\n")
cat("p-value:", binom_dfr$p.value, "\n")
if (binom_dfr$p.value < 0.05) {
  cat("Result: SIGNIFICANTLY HIGHER than expected\n\n")
} else {
  cat("Result: Not significantly higher\n\n")
}

# ============================================
# 4. Results
# ============================================

fisher_summary_v2 <- data.frame(
  Test = c("sul vs dfrA (Fisher)",
           "sul enrichment (Binomial)",
           "dfrA enrichment (Binomial)"),
  Observed = c(NA,
               round(mean(pa_matrix$sul_any) * 100, 1),
               round(mean(pa_matrix$dfrA_any) * 100, 1)),
  Expected = c(NA, expected_sul * 100, expected_dfrA * 100),
  p.value = c(fisher_sul_dfr$p.value,
              binom_sul$p.value,
              binom_dfr$p.value),
  Significant = c(fisher_sul_dfr$p.value < 0.05,
                  binom_sul$p.value < 0.05,
                  binom_dfr$p.value < 0.05),
  stringsAsFactors = FALSE
)

write.csv(fisher_summary_v2, 
          "~/amr_pilot/results/statistical_tests_v2.csv",
          row.names = FALSE)

cat("=== Analysis v2 Complete ===\n")
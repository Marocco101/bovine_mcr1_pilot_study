# ============================================
# Co-occurrence Analysis v3
# Author: Makiko FUJITA-SUZANNE
# Date: 10 April 2026
# Note: All 127 isolates are mcr-1 positive
# ============================================

# --- Read data ---
combined <- read.csv("~/amr_pilot/results/combined_resfinder_results.csv",
                     stringsAsFactors = FALSE)

str(combined)

# --- Define target genes ---
targets <- c("mcr-1.1", "sul1", "sul2", "dfrA1", "dfrA12", "dfrA17")

all_samples <- unique(combined$Sample)
cat("Number of isolates:", length(all_samples), "\n")

# ============================================
# Presence/Absence Matrix
# ============================================

pa_matrix <- data.frame(Sample = all_samples, stringsAsFactors = FALSE)

for (i in seq_along(targets)) {
  gene <- targets[i]
  col_name <- make.names(gene)
  samples_with_gene <- combined$Sample[combined$Resistance.gene == gene]
  pa_matrix[[gene]] <- as.integer(all_samples %in% samples_with_gene)
}

# --- mcr-1 all variants ---
all_genes <- unique(combined$Resistance.gene)
mcr1_variants <- all_genes[grep("^mcr-1", all_genes)]
cat("mcr-1 variants found:", mcr1_variants, "\n")

mcr1_samples <- combined$Sample[grep("^mcr-1", combined$Resistance.gene)]
pa_matrix$mcr1_any <- as.integer(all_samples %in% mcr1_samples)

# --- Aggregated columns ---
pa_matrix$sul_any <- as.integer(pa_matrix$sul1 == 1 | pa_matrix$sul2 == 1)

pa_matrix$dfrA_any <- as.integer(
  pa_matrix$dfrA1 == 1 | pa_matrix$dfrA12 == 1 | pa_matrix$dfrA17 == 1
)

pa_matrix$trap <- as.integer(pa_matrix$sul_any == 1 & pa_matrix$dfrA_any == 1)

# ============================================
# Detection Rates
# ============================================

cat("\n=== Dataset Summary ===\n")
cat("PRJNA1166088: mcr-1 positive E. coli from French cattle\n")
cat("Source: Haenni et al., 2025\n")
cat("N =", nrow(pa_matrix), "isolates\n")

cat("\n=== Detection Rates ===\n\n")

report_cols <- c("mcr-1.1", "sul1", "sul2", "sul_any",
                 "dfrA1", "dfrA12", "dfrA17", "dfrA_any", "trap")

detection <- sapply(report_cols, function(col) {
  n <- sum(pa_matrix[[col]])
  pct <- round(n / nrow(pa_matrix) * 100, 1)
  c(count = n, total = nrow(pa_matrix), percent = pct)
})

print(detection)

# ============================================
# Fisher's Exact Test: sul vs dfrA
# ============================================

cat("\n=== Fisher's Exact Test ===\n")
cat("Q: Do sul and dfrA co-occur more than expected by chance?\n\n")

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
  cat("variation for meaningful contingency table analysis.\n")
}

# ============================================
# Save results
# ============================================

write.csv(pa_matrix, "~/amr_pilot/results/presence_absence_matrix.csv",
          row.names = FALSE)

# --- Detection rates as CSV ---
det_df <- data.frame(
  Gene = report_cols,
  Count = as.integer(detection["count", ]),
  Total = as.integer(detection["total", ]),
  Percent = detection["percent", ],
  stringsAsFactors = FALSE
)

write.csv(det_df, "~/amr_pilot/results/detection_rates.csv",
          row.names = FALSE)

cat("\n=== Analysis Complete ===\n")
cat("Files saved:\n")
cat("  - presence_absence_matrix.csv\n")
cat("  - detection_rates.csv\n")
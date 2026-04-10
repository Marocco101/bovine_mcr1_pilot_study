# ============================================
# Co-occurrence Analysis & Fisher's Exact Test
# Author: Makiko FUJITA-SUZANNE
# Date: 9 April 2026
# ============================================

combined <- read.csv("~/amr_pilot/results/combined_resfinder_results.csv",
                     stringsAsFactors = FALSE)

str(combined)

targets <- c("mcr-1.1", "sul1", "sul2", "dfrA1", "dfrA12", "dfrA17")

all_samples <- unique(combined$Sample)
cat("Number of isolates:", length(all_samples), "\n")

# ============================================
# Presence/Absence Matrix
# ============================================

pa_matrix <- data.frame(Sample = all_samples, stringsAsFactors = FALSE)

for (i in seq_along(targets)) {
  gene <- targets[i]
  samples_with_gene <- combined$Sample[combined$Resistance.gene == gene]
  pa_matrix[[gene]] <- as.integer(all_samples %in% samples_with_gene)
}

all_genes <- unique(combined$Resistance.gene)
mcr1_variants <- all_genes[grep("^mcr-1", all_genes)]
cat("mcr-1 variants found:", mcr1_variants, "\n")

mcr1_samples <- combined$Sample[grep("^mcr-1", combined$Resistance.gene)]
pa_matrix$mcr1_any <- as.integer(all_samples %in% mcr1_samples)

pa_matrix$sul_any <- as.integer(pa_matrix$sul1 == 1 | pa_matrix$sul2 == 1)

pa_matrix$dfrA_any <- as.integer(
  pa_matrix$dfrA1 == 1 | pa_matrix$dfrA12 == 1 | pa_matrix$dfrA17 == 1
)

# ============================================
# Calculate Detection
# ============================================

cat("\n=== Detection Rates ===\n")

target_cols <- c("mcr1_any", "sul_any", "dfrA_any")

detection <- sapply(target_cols, function(col) {
  n <- sum(pa_matrix[[col]])
  pct <- round(n / nrow(pa_matrix) * 100, 1)
  c(count = n, total = nrow(pa_matrix), percent = pct)
})

print(detection)

# ============================================
# Fisher's Exact Test
# ============================================

run_fisher <- function(var1, var2, label) {
  tbl <- table(
    factor(var1, levels = c(0, 1)),
    factor(var2, levels = c(0, 1))
  )
  result <- fisher.test(tbl)
  
  cat("\n---", label, "---\n")
  print(tbl)
  cat("p-value:", result$p.value, "\n")
  cat("Odds ratio:", result$estimate, "\n")
  
  if (result$p.value < 0.05) {
    cat("Result: SIGNIFICANT (p < 0.05)\n")
  } else {
    cat("Result: Not significant\n")
  }
  
  return(list(table = tbl, p.value = result$p.value, 
              odds.ratio = result$estimate))
}

cat("\n=== Fisher's Exact Tests ===\n")

# Test 1: mcr-1 vs sul
result_sul <- run_fisher(pa_matrix$mcr1_any, pa_matrix$sul_any,
                         "mcr-1 vs sul")

# Test 2: mcr-1 vs dfrA
result_dfr <- run_fisher(pa_matrix$mcr1_any, pa_matrix$dfrA_any,
                         "mcr-1 vs dfrA")

# Test 3: mcr-1 vs sul&dfrA
pa_matrix$trap <- as.integer(pa_matrix$sul_any == 1 & pa_matrix$dfrA_any == 1)

result_trap <- run_fisher(pa_matrix$mcr1_any, pa_matrix$trap,
                          "mcr-1 vs sul&dfrA")

# ============================================
# Results Summary
# ============================================

fisher_summary <- data.frame(
  Test = c("mcr-1 vs sul", "mcr-1 vs dfrA", "mcr-1 vs Genetic Trap"),
  p.value = c(result_sul$p.value, result_dfr$p.value, result_trap$p.value),
  odds.ratio = c(result_sul$odds.ratio, result_dfr$odds.ratio, 
                 result_trap$odds.ratio),
  stringsAsFactors = FALSE
)

write.csv(pa_matrix, "~/amr_pilot/results/presence_absence_matrix.csv",
          row.names = FALSE)
write.csv(fisher_summary, "~/amr_pilot/results/fisher_test_results.csv",
          row.names = FALSE)

cat("\n=== Analysis Complete ===\n")
cat("Files saved:\n")
cat("  - ~/amr_pilot/results/presence_absence_matrix.csv\n")
cat("  - ~/amr_pilot/results/fisher_test_results.csv\n")
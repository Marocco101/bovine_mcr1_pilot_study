# ============================================
# ResFinder Results Integration
# Author: Makiko FUJITA-SUZANNE
# Date: 7 April 2026
# ============================================

results_dir <- "~/amr_pilot/results"

folders <- list.dirs(results_dir, recursive = FALSE)

all_results <- list()

for (i in seq_along(folders)) {
  sample <- basename(folders[i])
  file <- file.path(folders[i], "ResFinder_results_tab.txt")
  
  if (file.exists(file)) {
    df <- read.delim(file, sep = "\t", stringsAsFactors = FALSE)
    if (nrow(df) > 0) {
      df$Sample <- sample
      all_results[[sample]] <- df
    }
  }
}

combined <- do.call(rbind, all_results)
rownames(combined) <- NULL

cat("Total isolates processed:", length(all_results), "\n")
cat("Total gene hits:", nrow(combined), "\n")
cat("Unique resistance genes:", length(unique(combined$Resistance.gene)), "\n")

head(combined[, c("Sample", "Resistance.gene", "Identity", "Phenotype")])

write.csv(combined, "~/amr_pilot/results/combined_resfinder_results.csv", 
          row.names = FALSE)
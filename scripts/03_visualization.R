# ============================================
# Visualization for Pilot Study (v4 results)
# Author: Makiko FUJITA-SUZANNE
# Date: 17 April 2026
# Baseline: RESAPATH 2011-2018
# ============================================

pa_matrix <- read.csv(
  "~/Documents/bovine_mcr1_pilot_study/results/presence_absence_matrix.csv",
  stringsAsFactors = FALSE)

dir.create(
  "~/Documents/bovine_mcr1_pilot_study/results/figures",
  showWarnings = FALSE)

out_dir <- "~/Documents/bovine_mcr1_pilot_study/results/figures"

col_mcr  <- "#8B4513"
col_sul  <- "#4A708B"
col_dfr  <- "#556B2F"
col_trap <- "#B8860B"

# ============================================
# Figure 1: Individual Gene Detection Rates
# ============================================

gene_cols <- c("mcr.1.1", "sul1", "sul2", "dfrA1", "dfrA12", "dfrA17")
gene_labels <- c("mcr-1.1", "sul1", "sul2", "dfrA1", "dfrA12", "dfrA17")

rates <- sapply(gene_cols, function(col) {
  round(sum(pa_matrix[[col]]) / nrow(pa_matrix) * 100, 1)
})
names(rates) <- gene_labels

bar_colors <- c(col_mcr, col_sul, col_sul, col_dfr, col_dfr, col_dfr)

png(file.path(out_dir, "fig1_detection_rates.png"),
    width = 800, height = 500, res = 150)

par(mar = c(7, 5, 4, 2))

bp <- barplot(rates,
              col = bar_colors,
              ylim = c(0, 115),
              ylab = "Detection Rate (%)",
              main = "AMR Gene Detection in mcr-1 Positive E. coli\n(PRJNA1166088, n=127)",
              las = 2,
              cex.names = 0.9,
              cex.lab = 1.0,
              border = "grey40")

text(bp, rates + 3, paste0(rates, "%"), cex = 0.7)

legend("topright",
       legend = c("Colistin resistance",
                  "Sulfonamide resistance",
                  "Trimethoprim resistance"),
       fill = c(col_mcr, col_sul, col_dfr),
       border = "grey40",
       cex = 0.7,
       bg = "white")

dev.off()
cat("Figure 1 saved\n")

# ============================================
# Figure 2: Co-occurrence with RESAPATH baseline
# ============================================

# RESAPATH baseline values (2011-2018 mean) 
baseline_sul    <- 80.4
baseline_tmp    <- 30.1
baseline_tmpsul <- 38.0

categories <- c("sul (any)", "dfrA (any)", "Genetic Trap\n(sul + dfrA)")

obs_rates <- c(
  round(sum(pa_matrix$sul_any) / nrow(pa_matrix) * 100, 1),
  round(sum(pa_matrix$dfrA_any) / nrow(pa_matrix) * 100, 1),
  round(sum(pa_matrix$trap) / nrow(pa_matrix) * 100, 1)
)

baselines <- c(baseline_sul, baseline_tmp, baseline_tmpsul)

# Side-by-side bar plot
plot_data <- rbind(obs_rates, baselines)
rownames(plot_data) <- c("mcr-1+ isolates", "RESAPATH baseline")

bar_cols <- c(col_sul, "grey70")

png(file.path(out_dir, "fig2_resapath_comparison.png"),
    width = 800, height = 550, res = 150)

par(mar = c(8, 5, 4, 2))

bp2 <- barplot(plot_data,
               beside = TRUE,
               col = bar_cols,
               ylim = c(0, 115),
               ylab = "Prevalence (%)",
               main = "mcr-1+ Isolates vs RESAPATH Baseline\n(Bovine E. coli, 2011-2018)",
               names.arg = categories,
               las = 2,
               cex.names = 0.8,
               cex.lab = 1.0,
               border = "grey40")

# Add values on top of bars
text(bp2[1, ], obs_rates + 3,
     paste0(obs_rates, "%"), cex = 0.6, font = 2)
text(bp2[2, ], baselines + 3,
     paste0(baselines, "%"), cex = 0.6)

# Add significance stars
text(bp2[1, ], obs_rates + 9, "***", cex = 0.8, font = 2)

legend("topright",
       legend = c("mcr-1+ isolates (n=127)",
                  "RESAPATH bovine E. coli (2011-2018)"),
       fill = bar_cols,
       border = "grey40",
       cex = 0.65,
       bg = "white")

dev.off()
cat("Figure 2 saved\n")

# ============================================
# Figure 3: Heatmap (Presence/Absence Matrix)
# ============================================

heat_cols <- c("mcr.1.1", "sul1", "sul2", "dfrA1", "dfrA12", "dfrA17")
heat_labels <- c("mcr-1.1", "sul1", "sul2", "dfrA1", "dfrA12", "dfrA17")

heat_data <- as.matrix(pa_matrix[, heat_cols])
colnames(heat_data) <- heat_labels

# Sort by Trap status
trap_order <- order(pa_matrix$trap, pa_matrix$sul_any,
                    pa_matrix$dfrA_any, decreasing = TRUE)
heat_data <- heat_data[trap_order, ]

png(file.path(out_dir, "fig3_heatmap.png"),
    width = 800, height = 550, res = 150)

par(mar = c(6, 8, 4, 2))

image(t(heat_data[nrow(heat_data):1, ]),
      col = c("#F5F5F0", "#4A708B"),
      axes = FALSE,
      main = "Presence/Absence of Target AMR Genes\n(127 mcr-1+ isolates, sorted by Trap status)")

axis(1, at = seq(0, 1, length.out = ncol(heat_data)),
     labels = heat_labels, las = 2, cex.axis = 0.9)
axis(2, at = c(0, 1),
     labels = c(paste0("Isolate ", nrow(heat_data)), "Isolate 1"),
     las = 1, cex.axis = 0.7)

legend("bottomright",
       legend = c("Present", "Absent"),
       fill = c("#4A708B", "#F5F5F0"),
       border = "grey50",
       cex = 0.7,
       bg = "white")

dev.off()
cat("Figure 3 saved\n")

# ============================================
# Output
# ============================================
cat("\n=== All Figures Complete ===\n")
cat("Files:\n")
figs <- list.files(out_dir)
for (f in figs) {
  cat("  -", f, "\n")
}

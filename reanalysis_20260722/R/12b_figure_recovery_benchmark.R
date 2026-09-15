script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE)), "_common.R"))

sum <- read.csv(p("outputs", "recovery_benchmark", "recovery_benchmark_summary.csv"),
                check.names = FALSE)

out_dir <- p("outputs", "recovery_benchmark")
ensure_dirs(out_dir)

method_short <- c(
  naive_fixed_set_gene_label      = "Naive fixed set",
  gene_label_complete_rule_replay = "Gene-label replay",
  block_preserving_rule_replay    = "Block-preserving replay",
  full_two_arm_synthetic_replay   = "Full two-arm replay",
  sample_splitting_fixed_set      = "Sample splitting"
)
method_col <- c(
  naive_fixed_set_gene_label      = "#B2182B",
  gene_label_complete_rule_replay = "#2166AC",
  block_preserving_rule_replay    = "#4393C3",
  full_two_arm_synthetic_replay   = "#92C5DE",
  sample_splitting_fixed_set      = "#009E73"
)
conf_lty <- c(none = 1, adjusted = 2)

open_device <- function(device, file, width, height) {
  if (device == "pdf") pdf(file, width = width, height = height, useDingbats = FALSE)
  if (device == "png") png(file, width = width * 300, height = height * 300, res = 300)
  if (device == "tiff") tiff(file, width = width * 300, height = height * 300,
                             res = 300, compression = "lzw")
  par(fig = c(0, 1, 0, 1))
  plot.new()
}

panel_begin <- function(fig, main = NULL, ylab = NULL) {
  par(fig = fig, new = TRUE, mar = c(4.0, 4.2, 2.8, 0.8), las = 1,
      cex.axis = 0.85, cex.lab = 0.95, mgp = c(2.4, 0.7, 0))
}

panel_label <- function(label) {
  mtext(label, side = 3, adj = 0, line = 1.6, font = 2, cex = 1.0)
}

alpha_line <- function() {
  abline(h = 0.05, col = "#595959", lty = 2, lwd = 0.8)
}

devs <- c("pdf", "png", "tiff")
for (dev in devs) {
  open_device(dev, file.path(out_dir, paste0("Figure2_recovery_benchmark.", dev)),
              width = 10.8, height = 8.4)

  # ---- Panel A: type I error under global and single-arm nulls ----
  null_rows <- sum[sum$mechanism %in% c("global_null", "single_arm_null"), ]
  methods <- names(method_short)
  scen <- c("global_null", "single_arm_null")
  scen_lab <- c("Global null\n(no signal in either arm)", "Single-arm null\n(signal in one arm)")
  confs <- c("none", "adjusted")
  ncell <- length(methods) * length(confs)
  key_lookup <- interaction(null_rows$scenario, null_rows$method, null_rows$confounding)
  key_query <- interaction(rep(scen, each = ncell),
                           rep(rep(methods, each = length(confs)), times = length(scen)),
                           rep(confs, times = length(scen) * length(methods)))
  vals <- null_rows$rejection_rate[match(key_query, key_lookup)]
  lo <- null_rows$wilson_lower[match(key_query, key_lookup)]
  hi <- null_rows$wilson_upper[match(key_query, key_lookup)]
  adjusted_idx <- rep(confs, times = length(scen) * length(methods)) == "adjusted"

  panel_begin(c(0.04, 0.45, 0.53, 1.0), ylab = "Rejection rate (500 replicates)")
  bp <- barplot(matrix(vals, ncol = length(scen)), beside = TRUE,
                col = rep(method_col[methods], times = length(scen)),
                border = NA, ylim = c(0, 1.12), width = 0.85,
                names.arg = scen_lab, cex.names = 0.8, space = c(0.25, 0.8),
                axes = FALSE, xlab = NA,
                main = "A  Calibration: type I error under two null constructions",
                cex.main = 0.95)
  xmid <- as.vector(bp)
  axis(2, at = seq(0, 1, 0.25), labels = seq(0, 1, 0.25), cex.axis = 0.85)
  segments(xmid, lo, xmid, hi, col = "#333333", lwd = 0.6)
  points(xmid[adjusted_idx], vals[adjusted_idx], pch = 21, col = "#333333",
         bg = "white", cex = 0.85, lwd = 0.8)
  alpha_line()
  legend("topright",
         legend = c(method_short, "No confounder (bar)", "Adjusted confounder (open circle)"),
         fill = c(method_col, NA, NA), border = c(rep(NA, 5), NA, NA),
         pch = c(rep(NA, 5), 15, 21), col = c(rep(NA, 5), "#555555", "#333333"),
         pt.bg = c(rep(NA, 5), NA, "white"), pt.cex = 1.0,
         cex = 0.56, bty = "n")
  box()
  panel_label("A")

  # ---- Panel B: proportional shared effects, power versus h ----
  prop <- sum[sum$mechanism == "proportional_shared" &
                sum$method %in% c("gene_label_complete_rule_replay",
                                  "full_two_arm_synthetic_replay"), ]
  hs <- sort(unique(prop$h))
  panel_begin(c(0.52, 1.0, 0.53, 1.0), ylab = "Rejection rate")
  plot(NA, xlim = range(hs) + c(-0.15, 0.15), ylim = c(0, 1.05),
       xlab = "Shared-component magnitude h (within-arm SD units)", ylab = NA,
       xaxt = "n", yaxt = "n",
       main = "B  Recovery: rejection rate under proportional shared effects",
       cex.main = 0.95)
  axis(1, at = hs, labels = as.character(hs))
  axis(2, at = seq(0, 1, 0.25), labels = seq(0, 1, 0.25))
  for (cc in c(0.5, 1.0)) {
    for (conf in c("none", "adjusted")) {
      sub <- prop[prop$c == cc & prop$confounding == conf, ]
      for (m in c("gene_label_complete_rule_replay", "full_two_arm_synthetic_replay")) {
        row <- sub[sub$method == m, ]
        row <- row[order(row$h), ]
        lwd <- if (m == "gene_label_complete_rule_replay") 1.6 else 1.1
        col <- if (m == "gene_label_complete_rule_replay") "#2166AC" else "#92C5DE"
        lines(row$h, row$rejection_rate, type = "b",
              pch = if (m == "gene_label_complete_rule_replay") 19 else 4,
              lty = conf_lty[conf], col = col, lwd = lwd, cex = 0.9)
      }
    }
  }
  alpha_line()
  legend("topleft", legend = c("c = 0.5, gene replay", "c = 1.0, gene replay",
                               "c = 0.5, full two-arm", "c = 1.0, full two-arm"),
         col = rep(c("#2166AC", "#92C5DE"), each = 2), pch = c(19, 19, 4, 4),
         lty = rep(c(1, 2), 2), lwd = 1.2, cex = 0.6, bty = "n")
  legend("bottomright", legend = c("No confounder", "Adjusted confounder"),
         lty = c(1, 2), col = "#555555", lwd = 1.2, cex = 0.6, bty = "n")
  box()
  panel_label("B")

  # ---- Panel C: sparse and direction-heterogeneous mechanisms ----
  sparse <- sum[sum$mechanism == "sparse_shared" & sum$confounding == "none" &
                  sum$method %in% c("naive_fixed_set_gene_label",
                                    "gene_label_complete_rule_replay",
                                    "sample_splitting_fixed_set"), ]
  het <- sum[sum$mechanism == "heterogeneous" & sum$confounding == "none" &
               sum$method %in% c("naive_fixed_set_gene_label",
                                 "gene_label_complete_rule_replay",
                                 "block_preserving_rule_replay",
                                 "sample_splitting_fixed_set"), ]
  m1 <- matrix(sparse$rejection_rate, ncol = 2, byrow = FALSE)
  m2 <- matrix(het$rejection_rate, ncol = 2, byrow = FALSE)

  panel_begin(c(0.04, 0.34, 0.03, 0.47), ylab = "Rejection rate")
  barplot(m1, beside = TRUE, col = method_col[unique(sparse$method)],
          border = NA, ylim = c(0, 1.12), width = 0.7,
          names.arg = c("Sparse 0.8 SD\n(10% of genes)", "Sparse 1.2 SD\n(10% of genes)"),
          cex.names = 0.72, space = c(0.2, 0.5), axes = FALSE, xlab = NA,
          main = "C  Detection limits", cex.main = 0.95)
  axis(2, at = seq(0, 1, 0.25), labels = seq(0, 1, 0.25), cex.axis = 0.85)
  alpha_line()
  box()
  panel_label("C")

  panel_begin(c(0.38, 0.68, 0.03, 0.47), ylab = "Rejection rate")
  barplot(m2, beside = TRUE, col = method_col[unique(het$method)],
          border = NA, ylim = c(0, 1.12), width = 0.7,
          names.arg = c("25% flipped", "50% flipped"),
          cex.names = 0.75, space = c(0.2, 0.5), axes = FALSE, xlab = NA,
          main = "C (cont.)  Direction heterogeneity", cex.main = 0.95)
  axis(2, at = seq(0, 1, 0.25), labels = seq(0, 1, 0.25), cex.axis = 0.85)
  alpha_line()
  legend("topleft",
         legend = method_short[c("naive_fixed_set_gene_label",
                                 "gene_label_complete_rule_replay",
                                 "block_preserving_rule_replay",
                                 "sample_splitting_fixed_set")],
         fill = method_col[c("naive_fixed_set_gene_label",
                             "gene_label_complete_rule_replay",
                             "block_preserving_rule_replay",
                             "sample_splitting_fixed_set")],
         border = NA, cex = 0.56, bty = "n")
  box()
  panel_label("C")

  # ---- Panel D: conditional versus global null constructions ----
  key_scen <- c("global_null", "single_arm_null", "prop_shared_c0.5_h2.5",
                "prop_shared_c1_h1.5", "prop_shared_c1_h2.5")
  key_lab <- c("Global null", "Single-arm null", "c = 0.5, h = 2.5",
               "c = 1.0, h = 1.5", "c = 1.0, h = 2.5")
  d_rows <- sum[sum$scenario %in% key_scen & sum$confounding == "none" &
                  sum$method %in% c("gene_label_complete_rule_replay",
                                    "full_two_arm_synthetic_replay"), ]
  mat <- matrix(NA, nrow = 2, ncol = length(key_scen))
  rownames(mat) <- c("gene_label_complete_rule_replay", "full_two_arm_synthetic_replay")
  for (i in seq_along(key_scen)) {
    for (m in rownames(mat)) {
      row <- d_rows[d_rows$scenario == key_scen[i] & d_rows$method == m, ]
      if (nrow(row)) mat[m, i] <- row$rejection_rate
    }
  }
  panel_begin(c(0.72, 1.0, 0.03, 0.47), ylab = "Rejection rate")
  barplot(mat, beside = TRUE, col = c("#2166AC", "#92C5DE"), border = NA,
          ylim = c(0, 1.12), names.arg = key_lab, cex.names = 0.68,
          space = c(0.2, 0.5), axes = FALSE, xlab = NA,
          main = "D  Conditional one-arm vs global two-arm null", cex.main = 0.95)
  axis(2, at = seq(0, 1, 0.25), labels = seq(0, 1, 0.25), cex.axis = 0.85)
  alpha_line()
  legend("topleft", legend = c("Conditional one-arm (gene replay)",
                               "Global two-arm (full replay)"),
         fill = c("#2166AC", "#92C5DE"), border = NA, cex = 0.58, bty = "n")
  box()
  panel_label("D")

  dev.off()
}

message("Figure 2 (calibration-recovery benchmark) written to ",
        file.path(out_dir, "Figure2_recovery_benchmark.pdf"))

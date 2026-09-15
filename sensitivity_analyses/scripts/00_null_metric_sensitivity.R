# 00_null_metric_sensitivity.R
# Metric-sensitivity replay of the phenotype-level null (Freedman-Lane).
# Reuses the archived pipeline logic (reanalysis_pipeline/R/11_osa_phenotype_null.R)
# with the SAME seed prefix (50260723) so the permutation stream is reproducible
# in principle; a per-iteration comparison against the archived first 500
# iterations is reported, and distribution-level agreement (mean/sd/quantiles)
# is the validity check for the replayed null.
# Each iteration additionally computes Spearman and standardized-effect
# (log2FC / SE) correlations of the selected gene set. Diagnostic only.

script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
SCRIPT_DIR <- dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE))
BASE <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), winslash = "/")
# Current analysis tree; the reconstructed expression RDS is generated locally from public GEO inputs and is not redistributed.
ARCH <- file.path(BASE, "reanalysis_20260722")
OUT <- file.path(BASE, "sensitivity_analyses", "outputs")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages(library(limma))

seed <- 50260723L
n_perm <- 500L
fdr_cut <- 0.05
lfc_cut <- 0.3
min_selected <- 3L

y_chromosome_genes <- c(
  "UTY", "USP9Y", "DDX3Y", "KDM5D", "EIF1AY", "ZFY", "SRY",
  "NLGN4Y", "RPS4Y1", "RPS4Y2", "TSPY1", "RBMY1A1", "DAZ1",
  "PRKY", "AMELY", "TBL1Y", "PCDH11Y", "TMSB4Y", "VCY",
  "CDY1", "CDY2A", "HSFY1", "TXLNGY", "BPY2", "PRY"
)

safe_cor <- function(x, y, method = "pearson") {
  if (length(x) < min_selected || length(y) != length(x) ||
      stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  suppressWarnings(as.numeric(stats::cor(x, y, method = method)))
}

read_de <- function(path) {
  d <- read.csv(path, check.names = FALSE)
  out <- data.frame(
    Gene = trimws(as.character(d$Gene)),
    MASLD_log2FC = as.numeric(d$log2FoldChange),
    MASLD_FDR = as.numeric(d$padj),
    MASLD_SE = as.numeric(d$lfcSE),
    stringsAsFactors = FALSE
  )
  out[nzchar(out$Gene) & !duplicated(out$Gene), , drop = FALSE]
}

masld_de <- read_de(file.path(ARCH, "outputs", "de", "GSE126848_DESeq2_primary.csv"))
masld_de$MASLD_significant <- is.finite(masld_de$MASLD_FDR) &
  masld_de$MASLD_FDR < fdr_cut & abs(masld_de$MASLD_log2FC) > lfc_cut

masld_hubs <- unique(read.csv(
  file.path(ARCH, "outputs", "network", "MASLD_topology_hubs_primary.csv"),
  check.names = FALSE
)$Gene)
osa_hubs <- unique(read.csv(
  file.path(ARCH, "outputs", "osa_network", "OSA_topology_hubs_primary.csv"),
  check.names = FALSE
)$Gene)

osa_obj <- readRDS(file.path(ARCH, "data", "derived", "GSE135917_rebuilt.rds"))
osa_expr <- osa_obj$discovery_expression
osa_coldata <- osa_obj$discovery_coldata
stopifnot(identical(colnames(osa_expr), rownames(osa_coldata)), ncol(osa_expr) == 18L)

design_nuisance <- model.matrix(~ age_z + bmi_z + sex, data = osa_coldata)
design_full <- model.matrix(~ age_z + bmi_z + sex + group, data = osa_coldata)
fit_nuisance <- lmFit(osa_expr, design_nuisance)
nuisance_fitted <- fit_nuisance$coefficients %*% t(design_nuisance)
nuisance_residuals <- osa_expr - nuisance_fitted

fit_effect <- function(expression_matrix) {
  fit <- eBayes(lmFit(expression_matrix, design_full), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = "groupOSA", number = Inf, sort.by = "none")
  se <- sqrt(fit$s2.post) * fit$stdev.unscaled[, "groupOSA"]
  data.frame(
    Gene = rownames(tt),
    OSA_log2FC = as.numeric(tt$logFC),
    OSA_FDR = as.numeric(tt$adj.P.Val),
    OSA_SE = as.numeric(se),
    stringsAsFactors = FALSE
  )
}

evaluate_pipeline <- function(osa_de) {
  osa_de$OSA_significant <- is.finite(osa_de$OSA_FDR) &
    osa_de$OSA_FDR < fdr_cut & abs(osa_de$OSA_log2FC) > lfc_cut
  d <- merge(masld_de, osa_de, by = "Gene", all = FALSE, sort = FALSE)
  d <- d[is.finite(d$MASLD_log2FC) & is.finite(d$OSA_log2FC) &
           !d$Gene %in% y_chromosome_genes, , drop = FALSE]
  d$OSA_hub_and_MASLD_DE <- d$Gene %in% osa_hubs & d$MASLD_significant
  d$MASLD_hub_and_OSA_DE <- d$Gene %in% masld_hubs & d$OSA_significant
  pool <- d[d$OSA_hub_and_MASLD_DE | d$MASLD_hub_and_OSA_DE, , drop = FALSE]
  same <- pool$MASLD_log2FC * pool$OSA_log2FC > 0
  x <- pool$MASLD_log2FC[same]
  y <- pool$OSA_log2FC[same]
  r_pearson <- safe_cor(x, y, "pearson")
  r_spearman <- safe_cor(x, y, "spearman")
  xs <- x / pool$MASLD_SE[same]
  ys <- y / pool$OSA_SE[same]
  ok_std <- is.finite(xs) & is.finite(ys) & pool$OSA_SE[same] > 0 & pool$MASLD_SE[same] > 0
  r_std <- if (sum(ok_std) >= min_selected) {
    safe_cor(xs[ok_std], ys[ok_std], "pearson")
  } else {
    NA_real_
  }
  list(
    common = nrow(d),
    osa_significant = sum(osa_de$OSA_significant),
    eligible = nrow(pool),
    selected = sum(same),
    r_pearson = r_pearson,
    r_spearman = r_spearman,
    r_std = r_std
  )
}

observed <- evaluate_pipeline(fit_effect(osa_expr))
stopifnot(observed$common == 12958L, observed$eligible == 110L,
          observed$selected == 52L, is.finite(observed$r_pearson))
cat(sprintf("observed: n=%d Pearson=%.6f Spearman=%.6f Std=%.6f\n",
            observed$selected, observed$r_pearson, observed$r_spearman, observed$r_std))

set.seed(seed)
permutations <- replicate(n_perm, sample.int(ncol(osa_expr)), simplify = FALSE)
rows <- vector("list", n_perm)
start_time <- Sys.time()
for (i in seq_len(n_perm)) {
  rows[[i]] <- tryCatch({
    null_expr <- nuisance_fitted + nuisance_residuals[, permutations[[i]], drop = FALSE]
    result <- evaluate_pipeline(fit_effect(null_expr))
    data.frame(
      Iteration = i,
      Converged = TRUE,
      SameSignSelectedGenes = result$selected,
      PearsonR = result$r_pearson,
      SpearmanR = result$r_spearman,
      StdR = result$r_std,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      Iteration = i, Converged = FALSE, SameSignSelectedGenes = NA_integer_,
      PearsonR = NA_real_, SpearmanR = NA_real_, StdR = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  if (i %% 100L == 0L || i == n_perm) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    cat(sprintf("Metric-sensitivity null: %d/%d iterations (%.1f min elapsed)\n",
                i, n_perm, elapsed))
  }
  if (i %% 250L == 0L) gc(verbose = FALSE)
}

iterations <- do.call(rbind, rows)
write.csv(iterations, file.path(OUT, "null_metric_sensitivity_iterations.csv"),
          row.names = FALSE, na = "NA")

arch_it <- read.csv(file.path(ARCH, "outputs", "phenotype_null", "phenotype_null_iterations.csv"),
                    check.names = FALSE)
valid <- iterations$Converged & is.finite(iterations$PearsonR)
n_valid <- sum(valid)
arch_500 <- arch_it$SelectedPearson[seq_len(n_valid)]
arch_all <- arch_it$SelectedPearson
n_match <- sum(abs(iterations$PearsonR[valid] - arch_500) < 1e-12, na.rm = TRUE)
arch_500 <- arch_500[is.finite(arch_500)]
arch_all <- arch_all[is.finite(arch_all)]

summ_lines <- c(
  "NULL METRIC SENSITIVITY (500-iteration Freedman-Lane replay, seed prefix 50260723)",
  sprintf("valid iterations: %d; per-iteration Pearson agreement with archived first-%d iterations: %d/%d",
          n_valid, n_valid, n_match, n_valid),
  "  per-iteration identity is NOT expected to be bit-exact: the archived run's session state",
  "  (RNGkind / package versions / BLAS) is not archived (no session-info file for 11_osa_phenotype_null),",
  "  so the permutation stream cannot be reproduced bit-for-bit; validity is assessed at the",
  "  distribution level (below), which is the inference-relevant property.",
  sprintf("  distribution comparison, Pearson r under null (replay n=500 vs archived n=%d):",
          length(arch_all)),
  sprintf("    replay : mean=%.4f sd=%.4f median=%.4f q95=%.4f q975=%.4f",
          mean(iterations$PearsonR[valid], na.rm = TRUE),
          stats::sd(iterations$PearsonR[valid], na.rm = TRUE),
          stats::median(iterations$PearsonR[valid], na.rm = TRUE),
          unname(stats::quantile(iterations$PearsonR[valid], 0.95, na.rm = TRUE)),
          unname(stats::quantile(iterations$PearsonR[valid], 0.975, na.rm = TRUE))),
  sprintf("    archived: mean=%.4f sd=%.4f median=%.4f q95=%.4f q975=%.4f",
          mean(arch_all), stats::sd(arch_all), stats::median(arch_all),
          unname(stats::quantile(arch_all, 0.95)), unname(stats::quantile(arch_all, 0.975))),
  sprintf("    archived first %d: mean=%.4f sd=%.4f", n_valid, mean(arch_500), stats::sd(arch_500)),
  sprintf("observed: Pearson %.4f; Spearman %.4f; standardized-effect %.4f (n=%d)",
          observed$r_pearson, observed$r_spearman, observed$r_std, observed$selected),
  ""
)
obs_map <- c(PearsonR = observed$r_pearson, SpearmanR = observed$r_spearman, StdR = observed$r_std)
for (mcol in c("PearsonR", "SpearmanR", "StdR")) {
  v <- iterations[[mcol]][valid]
  v <- v[is.finite(v)]
  if (length(v) == 0) {
    summ_lines <- c(summ_lines, sprintf("%-12s no finite values", mcol))
    next
  }
  upper <- sum(v >= obs_map[[mcol]])
  p_up <- (upper + 1) / (length(v) + 1)
  summ_lines <- c(summ_lines, sprintf(
    "%-12s mean=%.4f sd=%.4f median=%.4f q95=%.4f P_upper(add-one)=%.4f",
    mcol, mean(v), stats::sd(v), stats::median(v),
    unname(stats::quantile(v, 0.95)), p_up
  ))
}

writeLines(summ_lines, file.path(OUT, "null_metric_sensitivity_summary.txt"))
cat(paste(summ_lines, collapse = "\n"), "\n")
cat("Metric-sensitivity null completed.\n")

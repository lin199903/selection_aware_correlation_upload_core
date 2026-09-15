# 13_osa_spike_in_recovery.R
# Spike-in recovery on the real GSE135917 residual covariance structure.
# Purpose: quantify the detection boundary of the conditional one-arm replay
# when shared effects are injected into the real transcriptome, rather than
# into the 120-gene Gaussian screen used by 12_osa_recovery_benchmark.R.
#
# Design
#   - Data: GSE135917 Study Group 1 discovery cohort (18 samples); nuisance
#     effects (age, BMI, sex) regressed out exactly as in 11_osa_phenotype_null.R.
#   - Injection set: n_spike genes drawn without replacement from the observed
#     topology-version eligible pool (110 genes); per outer replicate a fresh
#     draw, so the reported detection rate averages over shared-effect location.
#   - Injection: for each spiked gene, OSA-case samples get +delta/2 and
#     controls -delta/2 multiplied by sign(MASLD_log2FC); delta is in units of
#     one member gene's expression SD (z-score units).
#   - Null: Freedman-Lane residual permutation of the spiked matrix (nuisance
#     model refit on the spiked data; residual vectors permuted across samples,
#     preserving the gene covariance), pooled across all outer replicates.
#   - Outcome per dose: distribution of add-one P values across outer
#     replicates; rejection rate at alpha = 0.05.

script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE)), "_common.R"))

required_packages(c("limma", "jsonlite"))
suppressPackageStartupMessages(library(limma))

spec <- read_spec()
ensure_dirs(p("outputs", "spike_in_recovery"), p("logs"))

seed_base <- as.integer(spec$seed_benchmark) + 40000L
n_ref_perm <- as.integer(spec$benchmarks$recovery_benchmark$reference_replicates) * 2L
fdr_cut <- as.numeric(spec$de_thresholds$fdr)
lfc_cut <- as.numeric(spec$de_thresholds$abs_log2fc)
min_selected <- as.integer(spec$rule_replay$min_selected_genes)
alpha <- as.numeric(spec$benchmarks$alpha)

doses <- c(0.15, 0.30, 0.45, 0.60)
n_spike <- 25L
n_outer <- 30L
n_perm_per_outer <- as.integer(n_ref_perm / n_outer)

y_chromosome_genes <- c(
  "UTY", "USP9Y", "DDX3Y", "KDM5D", "EIF1AY", "ZFY", "SRY",
  "NLGN4Y", "RPS4Y1", "RPS4Y2", "TSPY1", "RBMY1A1", "DAZ1",
  "PRKY", "AMELY", "TBL1Y", "PCDH11Y", "TMSB4Y", "VCY",
  "CDY1", "CDY2A", "HSFY1", "TXLNGY", "BPY2", "PRY"
)

safe_cor <- function(x, y) {
  if (length(x) < min_selected || length(y) != length(x) ||
      stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  suppressWarnings(as.numeric(stats::cor(x, y, method = "pearson")))
}

read_de <- function(path, effect_name, fdr_name) {
  d <- read.csv(path, check.names = FALSE)
  stopifnot(all(c("Gene", "log2FoldChange", "padj") %in% names(d)))
  out <- data.frame(
    Gene = trimws(as.character(d$Gene)),
    effect = as.numeric(d$log2FoldChange),
    fdr = as.numeric(d$padj),
    stringsAsFactors = FALSE
  )
  names(out)[2:3] <- c(effect_name, fdr_name)
  out[nzchar(out$Gene) & !duplicated(out$Gene), , drop = FALSE]
}

masld_de <- read_de(
  p("outputs", "de", "GSE126848_DESeq2_primary.csv"),
  "MASLD_log2FC", "MASLD_FDR"
)
masld_de$MASLD_significant <- is.finite(masld_de$MASLD_FDR) &
  masld_de$MASLD_FDR < fdr_cut & abs(masld_de$MASLD_log2FC) > lfc_cut

masld_hubs <- unique(read.csv(
  p("outputs", "network", "MASLD_topology_hubs_primary.csv"),
  check.names = FALSE
)$Gene)
osa_hubs <- unique(read.csv(
  p("outputs", "osa_network", "OSA_topology_hubs_primary.csv"),
  check.names = FALSE
)$Gene)

osa_obj <- readRDS(p("data", "derived", "GSE135917_rebuilt.rds"))
osa_expr <- osa_obj$discovery_expression
osa_coldata <- osa_obj$discovery_coldata
stopifnot(identical(colnames(osa_expr), rownames(osa_coldata)), ncol(osa_expr) == 18L)

design_nuisance <- model.matrix(~ age_z + bmi_z + sex, data = osa_coldata)
design_full <- model.matrix(~ age_z + bmi_z + sex + group, data = osa_coldata)

fit_effect <- function(expression_matrix) {
  fit <- eBayes(lmFit(expression_matrix, design_full), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = "groupOSA", number = Inf, sort.by = "none")
  data.frame(
    Gene = rownames(tt),
    OSA_log2FC = as.numeric(tt$logFC),
    OSA_FDR = as.numeric(tt$adj.P.Val),
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
  list(
    common = nrow(d),
    osa_significant = sum(osa_de$OSA_significant),
    eligible = nrow(pool),
    selected = sum(same),
    r = safe_cor(pool$MASLD_log2FC[same], pool$OSA_log2FC[same])
  )
}

observed <- evaluate_pipeline(fit_effect(osa_expr))
stopifnot(observed$common == 12958L, observed$eligible == 110L,
          observed$selected == 52L, is.finite(observed$r))

pool_genes <- unique(read.csv(
  p("outputs", "primary", "topology_primary_eligible_pool.csv"),
  check.names = FALSE
)$Gene)
stopifnot(length(pool_genes) == 110L)

group_osa <- osa_coldata$group == "OSA"
stopifnot(sum(group_osa) > 0, sum(!group_osa) > 0)

masld_sign <- setNames(sign(masld_de$MASLD_log2FC), masld_de$Gene)

make_spiked <- function(expr, genes, delta, sign_vec) {
  out <- expr
  g <- genes[genes %in% rownames(out)]
  if (length(g) == 0L) return(out)
  s <- sign_vec[g]
  out[g, group_osa] <- out[g, group_osa, drop = FALSE] + delta / 2 * s
  out[g, !group_osa] <- out[g, !group_osa, drop = FALSE] - delta / 2 * s
  out
}

null_rows <- vector("list", n_outer)
obs_rows <- vector("list", n_outer)
gene_sets <- vector("list", n_outer)

for (j in seq_len(n_outer)) {
  set.seed(seed_base + j)
  spiked_genes <- sort(sample(pool_genes, n_spike, replace = FALSE))
  gene_sets[[j]] <- spiked_genes
  dose_results <- list()

  for (dose in doses) {
    spiked <- make_spiked(osa_expr, spiked_genes, dose, masld_sign)
    fit_n <- lmFit(spiked, design_nuisance)
    fitted_n <- fit_n$coefficients %*% t(design_nuisance)
    resid_n <- spiked - fitted_n
    obs <- evaluate_pipeline(fit_effect(spiked))
    obs_rows[[j]][[as.character(dose)]] <- obs

    set.seed(seed_base + 100000L + j)
    perms <- replicate(n_perm_per_outer, sample.int(ncol(spiked)), simplify = FALSE)
    rs <- vapply(perms, function(perm) {
      null_expr <- fitted_n + resid_n[, perm, drop = FALSE]
      evaluate_pipeline(fit_effect(null_expr))$r
    }, numeric(1))
    null_rows[[j]][[as.character(dose)]] <- rs
  }
  if (j %% 5L == 0L) cat(sprintf("Spike-in outer replicate %d/%d done\n", j, n_outer))
}

dose_summary <- lapply(doses, function(dose) {
  key <- as.character(dose)
  obs_r <- vapply(seq_len(n_outer), function(j) obs_rows[[j]][[key]]$r, numeric(1))
  null_all <- unlist(lapply(seq_len(n_outer), function(j) null_rows[[j]][[key]]))
  null_valid <- null_all[is.finite(null_all)]
  pvals <- vapply(obs_r, function(r) (1 + sum(null_valid >= r)) / (1 + length(null_valid)), numeric(1))
  data.frame(
    Dose = dose,
    NOuter = n_outer,
    NReferenceDraws = length(null_valid),
    MedianObservedR = median(obs_r, na.rm = TRUE),
    MinObservedR = min(obs_r, na.rm = TRUE),
    MaxObservedR = max(obs_r, na.rm = TRUE),
    MedianP = median(pvals, na.rm = TRUE),
    MinP = min(pvals, na.rm = TRUE),
    MaxP = max(pvals, na.rm = TRUE),
    RejectionRateAlpha05 = mean(pvals < alpha, na.rm = TRUE),
    NullMeanPearson = mean(null_valid, na.rm = TRUE),
    NullSDPearson = stats::sd(null_valid, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
summary_df <- do.call(rbind, dose_summary)
rownames(summary_df) <- NULL

write_csv_atomic(summary_df, p("outputs", "spike_in_recovery", "spike_in_recovery_summary.csv"))
saveRDS(list(
  doses = doses, n_spike = n_spike, n_outer = n_outer,
  spiked_gene_sets = gene_sets,
  observed = obs_rows, null = null_rows, summary = summary_df
), p("outputs", "spike_in_recovery", "spike_in_recovery_distribution.rds"), compress = "xz")

plot_detection <- function(device, file) {
  if (device == "pdf") pdf(file, width = 6.5, height = 5)
  if (device == "png") png(file, width = 1950, height = 1500, res = 300)
  on.exit(dev.off(), add = TRUE)
  plot(summary_df$Dose, summary_df$RejectionRateAlpha05, type = "b", pch = 19,
       col = "#2166AC", xlab = expression("Injected shared effect " * delta *
         " (member-gene z-score units)"), ylab = "Replay rejection rate at alpha = 0.05",
       ylim = c(0, 1), xlim = range(doses))
  abline(h = alpha, lty = 2, col = "grey40")
  text(max(doses) * 0.72, alpha + 0.05, labels = "nominal alpha", col = "grey40", cex = 0.9)
}
plot_detection("pdf", p("outputs", "spike_in_recovery", "spike_in_recovery_detection.pdf"))
plot_detection("png", p("outputs", "spike_in_recovery", "spike_in_recovery_detection.png"))

contract <- list(
  created = format(Sys.time(), tz = "UTC", usetz = TRUE),
  spec_version = spec$spec_version,
  accession = "GSE135917",
  population = "Study Group 1 case-control cohort",
  covariance = "observed residual covariance of the GSE135917 discovery matrix (age, BMI, sex regressed out)",
  injection = sprintf("%d genes drawn without replacement from the observed topology eligible pool (%d genes), sign-aligned with MASLD log2FC", n_spike, length(pool_genes)),
  effect_units = "one member gene expression SD (z-score units); case +delta/2, control -delta/2",
  null = "Freedman-Lane residual permutation of the spiked matrix (nuisance refit on spiked data; residual sample vectors permuted)",
  reference_draws_per_dose = n_perm_per_outer * n_outer,
  outer_replicates = n_outer,
  alpha = alpha,
  tail = "upper",
  add_one_correction = TRUE,
  primary = as.list(summary_df)
)
write_json_atomic(contract, p("outputs", "spike_in_recovery", "spike_in_recovery_contract.json"))

save_session_info("13_osa_spike_in_recovery")
print(summary_df)
cat("Spike-in recovery on real covariance completed.\n")

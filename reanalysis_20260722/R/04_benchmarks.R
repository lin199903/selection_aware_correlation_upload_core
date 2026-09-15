script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE)), "_common.R"))

required_packages(c("DESeq2", "limma", "jsonlite", "digest"))
suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
})

spec <- read_spec()
ensure_dirs(p("outputs", "benchmarks"), p("logs"))

seed <- as.integer(spec$seed_benchmark)
n_positive_perm <- as.integer(spec$rule_replay$n_permutations)
n_outer <- as.integer(spec$benchmarks$full_pipeline_null_iterations)
n_inner <- as.integer(spec$benchmarks$full_pipeline_null_inner_permutations)
alpha <- as.numeric(spec$benchmarks$alpha)
min_selected <- as.integer(spec$rule_replay$min_selected_genes)
fdr_cut <- as.numeric(spec$de_thresholds$fdr)
lfc_cut <- as.numeric(spec$de_thresholds$abs_log2fc)

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

add_one_p <- function(null, observed) {
  valid <- is.finite(null)
  if (!is.finite(observed) || !any(valid)) return(NA_real_)
  (1 + sum(null[valid] >= observed)) / (1 + sum(valid))
}

rule_replay <- function(x, y, n_perm, seed) {
  keep_obs <- x * y > 0
  observed_n <- sum(keep_obs)
  observed_r <- safe_cor(x[keep_obs], y[keep_obs])
  if (length(x) < min_selected || observed_n < min_selected || !is.finite(observed_r)) {
    return(list(observed_n = observed_n, observed_r = observed_r,
                null = numeric(), null_n = integer(), p_upper = NA_real_, valid = 0L))
  }
  set.seed(seed)
  null <- rep(NA_real_, n_perm)
  null_n <- integer(n_perm)
  for (i in seq_len(n_perm)) {
    yp <- sample(y, length(y), replace = FALSE)
    keep <- x * yp > 0
    null_n[i] <- sum(keep)
    if (null_n[i] >= min_selected) null[i] <- safe_cor(x[keep], yp[keep])
  }
  list(observed_n = observed_n, observed_r = observed_r,
       null = null, null_n = null_n, p_upper = add_one_p(null, observed_r),
       valid = sum(is.finite(null)))
}

read_de <- function(path, prefix) {
  d <- read.csv(path, check.names = FALSE)
  stopifnot(all(c("Gene", "log2FoldChange", "padj") %in% names(d)))
  out <- data.frame(
    Gene = trimws(as.character(d$Gene)),
    effect = as.numeric(d$log2FoldChange),
    fdr = as.numeric(d$padj),
    stringsAsFactors = FALSE
  )
  names(out)[2:3] <- paste0(prefix, c("_log2FC", "_FDR"))
  out[nzchar(out$Gene) & !duplicated(out$Gene), , drop = FALSE]
}

de_discovery <- read_de(p("outputs", "de", "GSE126848_DESeq2_primary.csv"), "Discovery")
de_replication <- read_de(p("outputs", "de", "GSE130970_DESeq2_replication.csv"), "Replication")
topology_hubs <- unique(read.csv(
  p("outputs", "network", "MASLD_topology_hubs_primary.csv"),
  check.names = FALSE
)$Gene)

# Independent empirical positive benchmark. The candidate set is fixed using
# discovery-cohort topology only; replication effects never define eligibility.
positive_common <- merge(de_discovery, de_replication, by = "Gene", sort = FALSE)
positive_common <- positive_common[
  is.finite(positive_common$Discovery_log2FC) &
    is.finite(positive_common$Replication_log2FC) &
    !positive_common$Gene %in% y_chromosome_genes,
  , drop = FALSE
]
positive_pool <- positive_common[positive_common$Gene %in% topology_hubs, , drop = FALSE]
positive <- rule_replay(
  positive_pool$Discovery_log2FC,
  positive_pool$Replication_log2FC,
  n_positive_perm,
  seed + 1L
)
positive_pool$SameSignObserved <- positive_pool$Discovery_log2FC *
  positive_pool$Replication_log2FC > 0
write_csv_atomic(positive_pool, p("outputs", "benchmarks", "independent_positive_pool.csv"))
write_csv_atomic(positive_pool[positive_pool$SameSignObserved, , drop = FALSE],
                 p("outputs", "benchmarks", "independent_positive_same_sign_selected.csv"))

positive_summary <- data.frame(
  Benchmark = "independent_empirical_positive",
  Discovery = "GSE126848 male-overlap liver",
  Replication = "GSE130970 independent liver",
  EligibilityRule = "GSE126848 topology hubs fixed before using replication effects",
  CommonGenes = nrow(positive_common),
  EligiblePoolGenes = nrow(positive_pool),
  SameSignSelectedGenes = positive$observed_n,
  WholeTranscriptomePearson = safe_cor(positive_common$Discovery_log2FC,
                                        positive_common$Replication_log2FC),
  WholeTranscriptomeSignAgreement = mean(
    sign(positive_common$Discovery_log2FC) == sign(positive_common$Replication_log2FC)
  ),
  ObservedSameSignPearson = positive$observed_r,
  RuleReplayNullMean = mean(positive$null, na.rm = TRUE),
  RuleReplayNullSD = sd(positive$null, na.rm = TRUE),
  RuleReplayPUpperAddOne = positive$p_upper,
  ValidPermutations = positive$valid,
  PositiveDetectedAtAlpha = is.finite(positive$p_upper) && positive$p_upper < alpha,
  stringsAsFactors = FALSE
)
write_csv_atomic(positive_summary,
                 p("outputs", "benchmarks", "independent_positive_summary.csv"))

# Conditional full-pipeline null on the rebuilt OSA discovery cohort. A
# Freedman-Lane construction removes age, BMI and sex, permutes whole-sample
# residual vectors (preserving gene-gene covariance), restores nuisance fits,
# and refits the OSA group effect. The fixed topology-only hub set and fixed
# MASLD evidence are then passed through cross-qualification, same-sign
# selection and the inner rule replay in every outer iteration.
osa_obj <- readRDS(p("data", "derived", "GSE135917_rebuilt.rds"))
osa_expr <- osa_obj$discovery_expression
osa_coldata <- osa_obj$discovery_coldata
stopifnot(identical(colnames(osa_expr), rownames(osa_coldata)), ncol(osa_expr) == 18L)

design_nuisance <- model.matrix(~ age_z + bmi_z + sex, data = osa_coldata)
design_full <- model.matrix(~ age_z + bmi_z + sex + group, data = osa_coldata)
fit_nuisance <- lmFit(osa_expr, design_nuisance)
nuisance_fitted <- fit_nuisance$coefficients %*% t(design_nuisance)
nuisance_residuals <- osa_expr - nuisance_fitted

masld_de <- data.frame(
  Gene = de_discovery$Gene,
  MASLD_log2FC = de_discovery$Discovery_log2FC,
  MASLD_FDR = de_discovery$Discovery_FDR,
  stringsAsFactors = FALSE
)
masld_de$MASLD_significant <- is.finite(masld_de$MASLD_FDR) &
  masld_de$MASLD_FDR < fdr_cut & abs(masld_de$MASLD_log2FC) > lfc_cut
osa_hubs <- unique(read.csv(
  p("outputs", "osa_network", "OSA_topology_hubs_primary.csv"),
  check.names = FALSE
)$Gene)

fit_permuted_de <- function(permutation) {
  null_expr <- nuisance_fitted + nuisance_residuals[, permutation, drop = FALSE]
  fit <- eBayes(lmFit(null_expr, design_full), trend = TRUE, robust = TRUE)
  res <- topTable(fit, coef = "groupOSA", number = Inf, sort.by = "none")
  data.frame(
    Gene = rownames(res), OSA_log2FC = as.numeric(res$logFC),
    OSA_FDR = as.numeric(res$adj.P.Val),
    stringsAsFactors = FALSE
  )
}

build_cross_disease_pool <- function(permuted_osa_de) {
  permuted_osa_de$OSA_significant <- is.finite(permuted_osa_de$OSA_FDR) &
    permuted_osa_de$OSA_FDR < fdr_cut & abs(permuted_osa_de$OSA_log2FC) > lfc_cut
  d <- merge(masld_de, permuted_osa_de, by = "Gene", all = FALSE, sort = FALSE)
  d <- d[is.finite(d$MASLD_log2FC) & is.finite(d$OSA_log2FC) &
           !d$Gene %in% y_chromosome_genes, , drop = FALSE]
  d$MASLD_significant <- is.finite(d$MASLD_FDR) & d$MASLD_FDR < fdr_cut &
    abs(d$MASLD_log2FC) > lfc_cut
  d$OSA_hub_and_MASLD_DE <- d$Gene %in% osa_hubs & d$MASLD_significant
  d$MASLD_hub_and_OSA_DE <- d$Gene %in% topology_hubs & d$OSA_significant
  d[d$OSA_hub_and_MASLD_DE | d$MASLD_hub_and_OSA_DE, , drop = FALSE]
}

outer_rows <- vector("list", n_outer)
set.seed(seed + 10000L)
residual_permutations <- replicate(n_outer, sample.int(ncol(osa_expr)), simplify = FALSE)
start_time <- Sys.time()
for (i in seq_len(n_outer)) {
  row <- tryCatch({
    perm_de <- fit_permuted_de(residual_permutations[[i]])
    pool <- build_cross_disease_pool(perm_de)
    replay <- rule_replay(pool$MASLD_log2FC, pool$OSA_log2FC,
                          n_inner, seed + 20000L + i)
    data.frame(
      Iteration = i,
      Converged = TRUE,
      SignificantOSAGenes = sum(is.finite(perm_de$OSA_FDR) &
                                  perm_de$OSA_FDR < fdr_cut &
                                  abs(perm_de$OSA_log2FC) > lfc_cut),
      EligiblePoolGenes = nrow(pool),
      SameSignSelectedGenes = replay$observed_n,
      ObservedPearson = replay$observed_r,
      RuleReplayPUpperAddOne = replay$p_upper,
      InnerValidPermutations = replay$valid,
      Error = NA_character_,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      Iteration = i, Converged = FALSE, SignificantOSAGenes = NA_integer_,
      EligiblePoolGenes = NA_integer_, SameSignSelectedGenes = NA_integer_,
      ObservedPearson = NA_real_, RuleReplayPUpperAddOne = NA_real_,
      InnerValidPermutations = NA_integer_, Error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
  outer_rows[[i]] <- row
  if (i %% 10L == 0L || i == n_outer) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    cat(sprintf("Full-pipeline null: %d/%d iterations (%.1f min elapsed)\n",
                i, n_outer, elapsed))
  }
  if (i %% 25L == 0L) gc(verbose = FALSE)
}
outer <- do.call(rbind, outer_rows)
write_csv_atomic(outer, p("outputs", "benchmarks", "full_pipeline_null_iterations.csv"))

valid_outer <- outer$Converged & is.finite(outer$RuleReplayPUpperAddOne)
n_valid_outer <- sum(valid_outer)
n_reject <- sum(outer$RuleReplayPUpperAddOne[valid_outer] < alpha)
rejection_rate <- if (n_valid_outer) n_reject / n_valid_outer else NA_real_

wilson_interval <- function(x, n, conf = 0.95) {
  if (n == 0L) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - conf) / 2)
  phat <- x / n
  den <- 1 + z^2 / n
  center <- (phat + z^2 / (2 * n)) / den
  half <- z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2)) / den
  c(max(0, center - half), min(1, center + half))
}
ci <- wilson_interval(n_reject, n_valid_outer)
negative_summary <- data.frame(
  Benchmark = "OSA_Freedman_Lane_conditional_full_pipeline_null",
  OuterIterationsPlanned = n_outer,
  OuterIterationsValid = n_valid_outer,
  FailedIterations = n_outer - n_valid_outer,
  InnerPermutationsPerIteration = n_inner,
  Alpha = alpha,
  Rejections = n_reject,
  RejectionRate = rejection_rate,
  Wilson95Lower = ci[1],
  Wilson95Upper = ci[2],
  AlphaInsideWilson95 = is.finite(ci[1]) && ci[1] <= alpha && alpha <= ci[2],
  MeanEligiblePoolGenes = mean(outer$EligiblePoolGenes[valid_outer], na.rm = TRUE),
  MeanSameSignSelectedGenes = mean(outer$SameSignSelectedGenes[valid_outer], na.rm = TRUE),
  stringsAsFactors = FALSE
)
write_csv_atomic(negative_summary,
                 p("outputs", "benchmarks", "full_pipeline_null_summary.csv"))

saveRDS(list(
  positive = positive,
  positive_pool = positive_pool,
  full_pipeline_null = outer,
  positive_summary = positive_summary,
  negative_summary = negative_summary
), p("outputs", "benchmarks", "benchmark_distributions.rds"), compress = "xz")

plot_benchmarks <- function(device, file) {
  if (device == "pdf") pdf(file, width = 12, height = 4.5)
  if (device == "png") png(file, width = 3600, height = 1350, res = 300)
  old <- par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
  on.exit({par(old); dev.off()}, add = TRUE)
  hist(positive$null, breaks = 60, col = "#D9EAF4", border = "white",
       main = "Independent positive benchmark", xlab = "Rule-replay Pearson r")
  abline(v = positive$observed_r, col = "#B2182B", lwd = 2, lty = 2)
  hist(outer$RuleReplayPUpperAddOne[valid_outer], breaks = seq(0, 1, by = 0.05),
       col = "#E8D9C5", border = "white", main = "OSA residual-null calibration",
       xlab = "Selection-aware P value")
  abline(v = alpha, col = "#B2182B", lwd = 2, lty = 2)
  plot(jitter(outer$EligiblePoolGenes[valid_outer]),
       jitter(outer$SameSignSelectedGenes[valid_outer]),
       pch = 16, cex = 0.7, col = grDevices::adjustcolor("#2166AC", alpha.f = 0.6),
       xlab = "Eligible pool genes", ylab = "Same-sign selected genes",
       main = "Outer-null set sizes")
}
plot_benchmarks("pdf", p("outputs", "benchmarks", "benchmark_diagnostics.pdf"))
plot_benchmarks("png", p("outputs", "benchmarks", "benchmark_diagnostics.png"))

contract <- list(
  created = format(Sys.time(), tz = "UTC", usetz = TRUE),
  spec_version = spec$spec_version,
  independent_positive = as.list(positive_summary),
  full_pipeline_null = as.list(negative_summary),
  full_pipeline_null_scope = "GSE135917 whole-sample nuisance residuals permuted; limma group effects, multiplicity threshold, cross-qualification, same-sign selection and inner rule replay repeated.",
  nuisance_estimation = "Freedman-Lane residual permutation under age_z + bmi_z + sex",
  nuisance_estimation_detail = "Nuisance fits are estimated once; one shared sample permutation is applied across all genes per outer iteration to preserve transcriptomic covariance.",
  network_step_recomputed_each_outer_iteration = FALSE,
  network_reuse_is_exact = FALSE,
  network_reuse_reason = "The benchmark is explicitly conditional on the rebuilt topology-only hub set; network estimation uncertainty is not included.",
  osa_side_status = "Rebuilt from GSE135917 Study Group 1 on GPL6244.",
  add_one_correction_used_for_all_inner_tests = TRUE
)

stopifnot(
  positive$valid > 0.99 * n_positive_perm,
  n_valid_outer >= 0.95 * n_outer,
  all(outer$RuleReplayPUpperAddOne[valid_outer] > 0 &
        outer$RuleReplayPUpperAddOne[valid_outer] <= 1)
)
write_json_atomic(contract, p("outputs", "benchmarks", "benchmark_contract.json"))
save_session_info("04_benchmarks")

print(positive_summary)
print(negative_summary)
cat("Independent positive and sample-label full-pipeline null benchmarks completed.\n")

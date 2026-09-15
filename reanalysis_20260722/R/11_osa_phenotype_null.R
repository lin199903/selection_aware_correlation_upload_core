script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE)), "_common.R"))

required_packages(c("limma", "jsonlite"))
suppressPackageStartupMessages(library(limma))

spec <- read_spec()
ensure_dirs(p("outputs", "phenotype_null"), p("logs"))

seed <- as.integer(spec$seed_benchmark) + 30000L
n_perm <- as.integer(spec$benchmarks$phenotype_null_iterations)
fdr_cut <- as.numeric(spec$de_thresholds$fdr)
lfc_cut <- as.numeric(spec$de_thresholds$abs_log2fc)
min_selected <- as.integer(spec$rule_replay$min_selected_genes)

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
fit_nuisance <- lmFit(osa_expr, design_nuisance)
nuisance_fitted <- fit_nuisance$coefficients %*% t(design_nuisance)
nuisance_residuals <- osa_expr - nuisance_fitted

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
      SignificantOSAGenes = result$osa_significant,
      CommonGenes = result$common,
      EligiblePoolGenes = result$eligible,
      SameSignSelectedGenes = result$selected,
      SelectedPearson = result$r,
      Error = NA_character_,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      Iteration = i, Converged = FALSE, SignificantOSAGenes = NA_integer_,
      CommonGenes = NA_integer_, EligiblePoolGenes = NA_integer_,
      SameSignSelectedGenes = NA_integer_, SelectedPearson = NA_real_,
      Error = conditionMessage(e), stringsAsFactors = FALSE
    )
  })
  if (i %% 100L == 0L || i == n_perm) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    cat(sprintf("Phenotype null: %d/%d iterations (%.1f min elapsed)\n",
                i, n_perm, elapsed))
  }
  if (i %% 250L == 0L) gc(verbose = FALSE)
}

iterations <- do.call(rbind, rows)
valid <- iterations$Converged & is.finite(iterations$SelectedPearson)
n_valid <- sum(valid)
n_upper <- sum(iterations$SelectedPearson[valid] >= observed$r)
p_upper <- (1 + n_upper) / (1 + n_valid)
mcse <- sqrt(p_upper * (1 - p_upper) / (n_valid + 1))
null_quantiles <- stats::quantile(
  iterations$SelectedPearson[valid], c(0.025, 0.5, 0.975), na.rm = TRUE
)

summary <- data.frame(
  Analysis = "OSA_Freedman_Lane_full_selection_replay",
  PlannedIterations = n_perm,
  ValidIterations = n_valid,
  FailedIterations = n_perm - n_valid,
  ObservedCommonGenes = observed$common,
  ObservedSignificantOSAGenes = observed$osa_significant,
  ObservedEligiblePoolGenes = observed$eligible,
  ObservedSameSignSelectedGenes = observed$selected,
  ObservedSelectedPearson = observed$r,
  NullMeanPearson = mean(iterations$SelectedPearson[valid]),
  NullSDPearson = stats::sd(iterations$SelectedPearson[valid]),
  NullQ025 = unname(null_quantiles[1]),
  NullMedian = unname(null_quantiles[2]),
  NullQ975 = unname(null_quantiles[3]),
  UpperTailCount = n_upper,
  PUpperAddOne = p_upper,
  MonteCarloSE = mcse,
  stringsAsFactors = FALSE
)

write_csv_atomic(iterations, p("outputs", "phenotype_null", "phenotype_null_iterations.csv"))
write_csv_atomic(summary, p("outputs", "phenotype_null", "phenotype_null_summary.csv"))
saveRDS(list(observed = observed, iterations = iterations, summary = summary),
        p("outputs", "phenotype_null", "phenotype_null_distribution.rds"), compress = "xz")

plot_null <- function(device, file) {
  if (device == "pdf") pdf(file, width = 7.2, height = 5.2)
  if (device == "png") png(file, width = 2160, height = 1560, res = 300)
  on.exit(dev.off(), add = TRUE)
  hist(iterations$SelectedPearson[valid], breaks = 50, col = "#D9EAF4",
       border = "white", xlab = "Selected-set Pearson r",
       main = "Phenotype-level residual null")
  abline(v = observed$r, col = "#B2182B", lwd = 2.5, lty = 2)
  legend("topleft", legend = sprintf("Observed r = %.3f; P = %.3f", observed$r, p_upper),
         lty = 2, lwd = 2.5, col = "#B2182B", bty = "n")
}
plot_null("pdf", p("outputs", "phenotype_null", "phenotype_null_distribution.pdf"))
plot_null("png", p("outputs", "phenotype_null", "phenotype_null_distribution.png"))

contract <- list(
  created = format(Sys.time(), tz = "UTC", usetz = TRUE),
  spec_version = spec$spec_version,
  accession = "GSE135917",
  population = "Study Group 1 case-control cohort",
  null = "Freedman-Lane residual permutation conditional on age, BMI and sex",
  permutation_unit = "whole-sample residual vector shared across genes",
  covariance_preserved_within_each_permuted_expression_matrix = TRUE,
  osa_group_effect_refit_each_iteration = TRUE,
  multiple_testing_and_cross_disease_selection_replayed_each_iteration = TRUE,
  network_step_recomputed_each_iteration = FALSE,
  network_conditioning = "MASLD and OSA topology-only hub sets held fixed",
  tail = "upper",
  add_one_correction = TRUE,
  primary_inference = as.list(summary)
)
write_json_atomic(contract, p("outputs", "phenotype_null", "phenotype_null_contract.json"))

stopifnot(n_valid >= 0.99 * n_perm, p_upper > 0, p_upper <= 1)
save_session_info("11_osa_phenotype_null")
print(summary)
cat("Phenotype-level residual null completed.\n")

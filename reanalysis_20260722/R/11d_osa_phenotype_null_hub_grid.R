script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE)), "_common.R"))

required_packages(c("limma", "jsonlite"))
suppressPackageStartupMessages(library(limma))

spec <- read_spec()
ensure_dirs(p("outputs", "phenotype_null_hub_grid"), p("logs"))

n_grid_perm <- 500L
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

masld_modules <- read.csv(p("outputs", "network", "MASLD_module_assignments.csv"),
                          check.names = FALSE)
osa_modules <- read.csv(p("outputs", "osa_network", "OSA_module_assignments.csv"),
                        check.names = FALSE)

build_hubs <- function(module_table, kme_cut, max_per_module) {
  rows <- lapply(setdiff(unique(module_table$Module), "grey"), function(mod) {
    d <- module_table[module_table$Module == mod &
                        is.finite(module_table$kME) &
                        abs(module_table$kME) >= kme_cut, , drop = FALSE]
    d <- d[order(-abs(d$kME), -d$MAD), , drop = FALSE]
    head(d, as.integer(max_per_module))
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) out <- module_table[0, ]
  unique(out$Gene)
}

grid <- data.frame(
  Definition = paste0("G", 1:8),
  KmeCut = c(0.70, 0.75, 0.80, 0.85, 0.90, 0.80, 0.80, 0.70),
  MaxPerModule = c(20L, 20L, 20L, 20L, 20L, 10L, 40L, 40L),
  stringsAsFactors = FALSE
)
grid$IsPrimary <- grid$KmeCut == 0.80 & grid$MaxPerModule == 20L

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

evaluate_pipeline <- function(osa_de, masld_hubs, osa_hubs) {
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

observed_full <- evaluate_pipeline(
  fit_effect(osa_expr),
  unique(read.csv(p("outputs", "network", "MASLD_topology_hubs_primary.csv"),
                  check.names = FALSE)$Gene),
  unique(read.csv(p("outputs", "osa_network", "OSA_topology_hubs_primary.csv"),
                  check.names = FALSE)$Gene)
)
stopifnot(observed_full$common == 12958L, observed_full$eligible == 110L,
          observed_full$selected == 52L, is.finite(observed_full$r))

grid_summary <- vector("list", nrow(grid))
grid_iterations <- vector("list", nrow(grid))
for (j in seq_len(nrow(grid))) {
  g <- grid[j, ]
  masld_hubs <- build_hubs(masld_modules, g$KmeCut, g$MaxPerModule)
  osa_hubs <- build_hubs(osa_modules, g$KmeCut, g$MaxPerModule)
  observed <- evaluate_pipeline(fit_effect(osa_expr), masld_hubs, osa_hubs)
  if (g$IsPrimary) {
    stopifnot(observed$eligible == 110L, observed$selected == 52L,
              is.finite(observed$r), length(masld_hubs) == 216L,
              length(osa_hubs) == 270L)
  }
  seed <- as.integer(spec$seed_benchmark) + 33000L + (j - 1L) * 1000L
  set.seed(seed)
  permutations <- replicate(n_grid_perm, sample.int(ncol(osa_expr)), simplify = FALSE)
  rows <- vector("list", n_grid_perm)
  start_time <- Sys.time()
  for (i in seq_len(n_grid_perm)) {
    null_expr <- nuisance_fitted + nuisance_residuals[, permutations[[i]], drop = FALSE]
    rows[[i]] <- tryCatch({
      result <- evaluate_pipeline(fit_effect(null_expr), masld_hubs, osa_hubs)
      data.frame(
        Definition = g$Definition,
        Iteration = i,
        Converged = TRUE,
        EligiblePoolGenes = result$eligible,
        SameSignSelectedGenes = result$selected,
        SelectedPearson = result$r,
        Error = NA_character_,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(
        Definition = g$Definition, Iteration = i, Converged = FALSE,
        EligiblePoolGenes = NA_integer_, SameSignSelectedGenes = NA_integer_,
        SelectedPearson = NA_real_, Error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    })
    if (i %% 100L == 0L || i == n_grid_perm) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      cat(sprintf("Hub grid %s (%s): %d/%d iterations (%.1f min elapsed)\n",
                  g$Definition, g$IsPrimary, i, n_grid_perm, elapsed))
    }
  }
  iters <- do.call(rbind, rows)
  valid <- iters$Converged & is.finite(iters$SelectedPearson)
  n_valid <- sum(valid)
  n_upper <- sum(iters$SelectedPearson[valid] >= observed$r)
  p_upper <- (1 + n_upper) / (1 + n_valid)
  mcse <- sqrt(p_upper * (1 - p_upper) / (n_valid + 1))
  q <- stats::quantile(iters$SelectedPearson[valid], c(0.025, 0.5, 0.975), na.rm = TRUE)
  grid_summary[[j]] <- data.frame(
    Definition = g$Definition,
    KmeCut = g$KmeCut,
    MaxPerModule = g$MaxPerModule,
    IsPrimary = g$IsPrimary,
    MasldHubs = length(masld_hubs),
    OsaHubs = length(osa_hubs),
    ObservedEligiblePoolGenes = observed$eligible,
    ObservedSameSignSelectedGenes = observed$selected,
    ObservedSelectedPearson = observed$r,
    ValidIterations = n_valid,
    FailedIterations = n_grid_perm - n_valid,
    NullMeanPearson = mean(iters$SelectedPearson[valid]),
    NullSDPearson = stats::sd(iters$SelectedPearson[valid]),
    NullQ025 = unname(q[1]),
    NullMedian = unname(q[2]),
    NullQ975 = unname(q[3]),
    UpperTailCount = n_upper,
    PUpperAddOne = p_upper,
    MonteCarloSE = mcse,
    stringsAsFactors = FALSE
  )
  grid_iterations[[j]] <- iters
  gc(verbose = FALSE)
}

summary <- do.call(rbind, grid_summary)
iterations <- do.call(rbind, grid_iterations)

write_csv_atomic(summary, p("outputs", "phenotype_null_hub_grid", "hub_grid_summary.csv"))
write_csv_atomic(iterations, p("outputs", "phenotype_null_hub_grid", "hub_grid_iterations.csv"))
saveRDS(list(summary = summary, iterations = iterations),
        p("outputs", "phenotype_null_hub_grid", "hub_grid_distributions.rds"),
        compress = "xz")

plot_grid <- function(device, file) {
  if (device == "pdf") pdf(file, width = 8.2, height = 5.6)
  if (device == "png") png(file, width = 2460, height = 1680, res = 300)
  on.exit(dev.off(), add = TRUE)
  col_primary <- "#B2182B"
  col_alt <- "#1B7837"
  plot(summary$KmeCut, summary$ObservedEligiblePoolGenes, type = "p",
       pch = ifelse(summary$IsPrimary, 17, 19),
       cex = ifelse(summary$IsPrimary, 1.6, 1.1),
       col = ifelse(summary$IsPrimary, col_primary, col_alt),
       xlab = "|kME| hub threshold", ylab = "Eligible pool size (genes)",
       main = "Hub-definition grid: eligible pool and replay P")
  text(summary$KmeCut, summary$ObservedEligiblePoolGenes,
       labels = sprintf("P=%.2f", summary$PUpperAddOne), pos = 3, cex = 0.6)
  abline(h = 110, lty = 3, col = "grey50")
  legend("topright",
         legend = c("Primary definition (|kME| >= 0.8, <= 20/module)", "Alternative definition"),
         pch = c(17, 19), col = c(col_primary, col_alt), bty = "n")
}
plot_grid("pdf", p("outputs", "phenotype_null_hub_grid", "hub_grid_pool_and_p.pdf"))
plot_grid("png", p("outputs", "phenotype_null_hub_grid", "hub_grid_pool_and_p.png"))

contract <- list(
  created = format(Sys.time(), tz = "UTC", usetz = TRUE),
  spec_version = spec$spec_version,
  accession = "GSE135917",
  population = "Study Group 1 case-control cohort",
  null = "Freedman-Lane residual permutation conditional on age, BMI and sex",
  hub_definitions = paste0(
    "non-grey modules; |kME| >= threshold; sorted by -abs(kME), -MAD; ",
    "capped at max genes per module; network held fixed across iterations"
  ),
  grid_permutations_per_definition = n_grid_perm,
  grid_permutation_rationale = "down-sampled relative to the primary 2000-iteration null (scripts 11/11c); Monte Carlo SE reported per row",
  primary_definition = "G3 (|kME| >= 0.80, max 20 per module); replicates primary analysis 110 eligible / 52 selected / r = 0.5507",
  tail = "upper",
  add_one_correction = TRUE,
  grid_inference = as.list(summary)
)
write_json_atomic(contract, p("outputs", "phenotype_null_hub_grid", "hub_grid_contract.json"))

stopifnot(all(summary$ValidIterations >= 0.99 * n_grid_perm),
          all(summary$PUpperAddOne > 0 & summary$PUpperAddOne <= 1))
save_session_info("11d_osa_phenotype_null_hub_grid")
print(summary)
cat("Hub-definition grid sensitivity completed.\n")

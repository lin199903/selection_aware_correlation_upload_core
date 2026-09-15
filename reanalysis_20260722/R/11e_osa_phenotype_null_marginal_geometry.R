script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE)), "_common.R"))

required_packages(c("limma", "jsonlite"))
suppressPackageStartupMessages(library(limma))

## 11e — Marginal-geometry diagnostic for the study-specific replay.
## Reproduces the exact draw sequence of 11_osa_phenotype_null.R (same seed,
## same replicate() call) and additionally records the per-draw selected gene
## identities with their MASLD / OSA effect estimates, so the marginal effect
## geometry of observed vs replayed selected sets can be compared.
## Chunked: args <chunk_start> <chunk_size> (1-based draw index).

args <- commandArgs(trailingOnly = TRUE)
chunk_start <- as.integer(args[1]); chunk_size <- as.integer(args[2])
stopifnot(is.finite(chunk_start), is.finite(chunk_size), chunk_start >= 1L, chunk_size >= 1L)

spec <- read_spec()
ensure_dirs(p("outputs", "phenotype_null_geometry"), p("logs"))

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

## Same pipeline as 11_, but returns the full selected-set detail.
evaluate_pipeline_detail <- function(osa_de) {
  osa_de$OSA_significant <- is.finite(osa_de$OSA_FDR) &
    osa_de$OSA_FDR < fdr_cut & abs(osa_de$OSA_log2FC) > lfc_cut
  d <- merge(masld_de, osa_de, by = "Gene", all = FALSE, sort = FALSE)
  d <- d[is.finite(d$MASLD_log2FC) & is.finite(d$OSA_log2FC) &
           !d$Gene %in% y_chromosome_genes, , drop = FALSE]
  d$OSA_hub_and_MASLD_DE <- d$Gene %in% osa_hubs & d$MASLD_significant
  d$MASLD_hub_and_OSA_DE <- d$Gene %in% masld_hubs & d$OSA_significant
  pool <- d[d$OSA_hub_and_MASLD_DE | d$MASLD_hub_and_OSA_DE, , drop = FALSE]
  same <- pool$MASLD_log2FC * pool$OSA_log2FC > 0
  sel <- pool[same, , drop = FALSE]
  list(
    common = nrow(d),
    osa_significant = sum(osa_de$OSA_significant),
    eligible = nrow(pool),
    selected = nrow(sel),
    genes = sel$Gene,
    masld_lfc = sel$MASLD_log2FC,
    osa_lfc = sel$OSA_log2FC
  )
}

## Observed selected set (iteration 0)
obs <- evaluate_pipeline_detail(fit_effect(osa_expr))
stopifnot(obs$common == 12958L, obs$eligible == 110L, obs$selected == 52L)
obs_r <- suppressWarnings(as.numeric(cor(obs$masld_lfc, obs$osa_lfc)))
stopifnot(abs(obs_r - 0.550719414222321) < 1e-12)
obs_df <- data.frame(
  Iteration = 0L, Gene = obs$genes,
  MASLD_log2FC = obs$masld_lfc, OSA_log2FC = obs$osa_lfc,
  stringsAsFactors = FALSE
)
write_csv_atomic(obs_df, p("outputs", "phenotype_null_geometry", "observed_selected_set.csv"))

## Regenerate the identical permutation sequence and take this chunk.
set.seed(seed)
permutations <- replicate(n_perm, sample.int(ncol(osa_expr)), simplify = FALSE)
idx <- chunk_start:min(chunk_start + chunk_size - 1L, n_perm)

rows <- vector("list", length(idx))
t0 <- Sys.time()
for (k in seq_along(idx)) {
  i <- idx[k]
  null_expr <- nuisance_fitted + nuisance_residuals[, permutations[[i]], drop = FALSE]
  res <- evaluate_pipeline_detail(fit_effect(null_expr))
  rows[[k]] <- data.frame(
    Iteration = i, Gene = res$genes,
    MASLD_log2FC = res$masld_lfc, OSA_log2FC = res$osa_lfc,
    stringsAsFactors = FALSE
  )
}
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("Chunk %d-%d done in %.2f min (%.2f s/draw)\n",
            chunk_start, max(idx), elapsed, 60 * elapsed / length(idx)))

out <- do.call(rbind, rows)
chunk_file <- p("outputs", "phenotype_null_geometry",
                sprintf("replay_selected_genes_%04d_%04d.csv", chunk_start, max(idx)))
write_csv_atomic(out, chunk_file)

## Verification against the archived iteration table for this chunk.
archived <- read.csv(p("outputs", "phenotype_null", "phenotype_null_iterations.csv"),
                     check.names = FALSE)
for (k in seq_along(idx)) {
  i <- idx[k]
  sel_k <- rows[[k]]
  a <- archived[archived$Iteration == i, ]
  stopifnot(nrow(a) == 1L)
  stopifnot(nrow(unique(sel_k["Iteration"])) == 1L)
  n_sel <- sum(sel_k$Iteration == i)
  stopifnot(n_sel == a$SameSignSelectedGenes)
  r_k <- suppressWarnings(as.numeric(cor(sel_k$MASLD_log2FC, sel_k$OSA_log2FC)))
  stopifnot(abs(r_k - a$SelectedPearson) < 1e-12)
}
cat("Chunk verified against archived iterations (counts and r identical).\n")

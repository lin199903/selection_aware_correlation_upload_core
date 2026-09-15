script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE)), "_common.R"))

required_packages(c("limma", "jsonlite", "parallel"))
suppressPackageStartupMessages({library(limma); library(parallel)})

## 12 — Same-disease MASLD stress test (historical filename retained for provenance).
## Primary: anchor GSE126848 (fixed DESeq2 effects) -> replay GSE135251 (limma on vst).
## Usage: Rscript 12_masld_masld_positive_control.R obs
##        Rscript 12_masld_masld_positive_control.R chunk <start> <size> <workers>

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
stopifnot(mode %in% c("obs", "chunk"))

spec <- read_spec()
ensure_dirs(p("outputs", "positive_control_masld_masld"), p("logs"))

seed <- as.integer(spec$seed_benchmark) + 31000L   # 20290724, frozen for this test
n_perm <- as.integer(spec$benchmarks$phenotype_null_iterations)  # 2000
fdr_cut <- as.numeric(spec$de_thresholds$fdr)
lfc_cut <- as.numeric(spec$de_thresholds$abs_log2fc)
min_selected <- as.integer(spec$rule_replay$min_selected_genes)

y_chromosome_genes <- c(
  "UTY", "USP9Y", "DDX3Y", "KDM5D", "EIF1AY", "ZFY", "SRY",
  "NLGN4Y", "RPS4Y1", "RPS4Y2", "TSPY1", "RBMY1A1", "DAZ1",
  "PRKY", "AMELY", "TBL1Y", "PCDH11Y", "TMSB4Y", "VCY",
  "CDY1", "CDY2A", "HSFY1", "TXLNGY", "BPY2", "PRY"
)

## ---- anchor arm: GSE126848, frozen DESeq2 effects ----
anchor_de <- read.csv(p("outputs", "de", "GSE126848_DESeq2_primary.csv"), check.names = FALSE)
stopifnot(all(c("Gene", "log2FoldChange", "padj") %in% names(anchor_de)))
anchor <- data.frame(
  Gene = trimws(as.character(anchor_de$Gene)),
  A_log2FC = as.numeric(anchor_de$log2FoldChange),
  A_FDR = as.numeric(anchor_de$padj),
  stringsAsFactors = FALSE
)
anchor <- anchor[nzchar(anchor$Gene) & !duplicated(anchor$Gene), , drop = FALSE]
anchor$A_significant <- is.finite(anchor$A_FDR) & anchor$A_FDR < fdr_cut &
  is.finite(anchor$A_log2FC) & abs(anchor$A_log2FC) > lfc_cut

hubs <- unique(trimws(read.csv(
  p("outputs", "network", "MASLD_topology_hubs_primary.csv"), check.names = FALSE)$Gene))
hubs <- hubs[nzchar(hubs)]

## ---- test arm: GSE135251, limma on vst (frozen substitution, contract section 4) ----
g251 <- readRDS(normalizePath(p("data", "derived", "GSE135251_rebuilt.rds"), winslash = "/", mustWork = TRUE))
expr <- g251$vst
coldata <- g251$metadata
stopifnot(identical(colnames(expr), coldata$sample_id), ncol(expr) == 216L)
coldata$group <- factor(coldata$group, levels = c("Control", "MASLD"))

design_full <- model.matrix(~ group, data = coldata)
design_nuisance <- model.matrix(~ 1, data = coldata)
fit_nuisance <- lmFit(expr, design_nuisance)
nuisance_fitted <- fit_nuisance$coefficients %*% t(design_nuisance)
nuisance_residuals <- expr - nuisance_fitted

fit_effect <- function(m) {
  fit <- eBayes(lmFit(m, design_full), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = "groupMASLD", number = Inf, sort.by = "none")
  data.frame(Gene = rownames(tt), B_log2FC = as.numeric(tt$logFC),
             B_FDR = as.numeric(tt$adj.P.Val), stringsAsFactors = FALSE)
}

evaluate_pipeline <- function(b_de) {
  b_de$B_significant <- is.finite(b_de$B_FDR) & b_de$B_FDR < fdr_cut &
    is.finite(b_de$B_log2FC) & abs(b_de$B_log2FC) > lfc_cut
  d <- merge(anchor, b_de, by = "Gene", all = FALSE, sort = FALSE)
  d <- d[is.finite(d$A_log2FC) & is.finite(d$B_log2FC) &
           !d$Gene %in% y_chromosome_genes, , drop = FALSE]
  d$eligible <- d$Gene %in% hubs & (d$A_significant | d$B_significant)
  pool <- d[d$eligible, , drop = FALSE]
  same <- pool$A_log2FC * pool$B_log2FC > 0
  sel <- pool[same, , drop = FALSE]
  r <- if (nrow(sel) < min_selected || sd(sel$A_log2FC) == 0 || sd(sel$B_log2FC) == 0) NA_real_
       else suppressWarnings(as.numeric(cor(sel$A_log2FC, sel$B_log2FC)))
  list(common = nrow(d), b_significant = sum(b_de$B_significant),
       eligible = nrow(pool), selected = nrow(sel), r = r,
       genes = sel$Gene, a_lfc = sel$A_log2FC, b_lfc = sel$B_log2FC)
}

if (mode == "obs") {
  obs <- evaluate_pipeline(fit_effect(expr))
  cat("observed: common", obs$common, "| eligible", obs$eligible,
      "| selected", obs$selected, "| r", obs$r, "\n")
  obs_df <- data.frame(Gene = obs$genes, A_log2FC = obs$a_lfc, B_log2FC = obs$b_lfc,
                       stringsAsFactors = FALSE)
  write_csv_atomic(obs_df, p("outputs", "positive_control_masld_masld", "observed_selected_set.csv"))

  ## whole-transcriptome context (descriptive only, NOT a success criterion)
  d_all <- merge(anchor, fit_effect(expr), by = "Gene", all = FALSE, sort = FALSE)
  d_all <- d_all[is.finite(d_all$A_log2FC) & is.finite(d_all$B_log2FC) &
                   !d_all$Gene %in% y_chromosome_genes, , drop = FALSE]
  wt_r <- cor(d_all$A_log2FC, d_all$B_log2FC)

  ## naive zero-centered comparator: 100,000 randomly paired sets of observed size
  set.seed(seed + 1L)
  n_sel <- obs$selected
  naive <- replicate(100000L, {
    i <- sample.int(nrow(d_all), n_sel)
    cor(d_all$A_log2FC[i], d_all$B_log2FC[i])
  })
  naive_p <- (1 + sum(abs(naive) >= abs(obs$r))) / (1 + length(naive))
  obs_summary <- data.frame(
    Analysis = "MASLD_MASLD_positive_control_observed",
    CommonGenes = obs$common, WholeTranscriptomePearson = wt_r,
    AnchorSignificant = sum(anchor$A_significant), TestArmSignificant = obs$b_significant,
    EligiblePoolGenes = obs$eligible, SameSignSelectedGenes = obs$selected,
    ObservedSelectedPearson = obs$r, NaiveP_two_sided = naive_p,
    NaiveDraws = length(naive), stringsAsFactors = FALSE)
  write_csv_atomic(obs_summary, p("outputs", "positive_control_masld_masld", "observed_summary.csv"))
  print(obs_summary)
  cat("Observed analysis completed.\n")
  quit(save = "no")
}

## ---- chunk mode ----
chunk_start <- as.integer(args[2]); chunk_size <- as.integer(args[3]); nworkers <- as.integer(args[4])
stopifnot(is.finite(chunk_start), is.finite(chunk_size), chunk_start >= 1L)

set.seed(seed)
permutations <- replicate(n_perm, sample.int(ncol(expr)), simplify = FALSE)
idx <- chunk_start:min(chunk_start + chunk_size - 1L, n_perm)

worker <- function(i) {
  null_expr <- nuisance_fitted + nuisance_residuals[, permutations[[i]], drop = FALSE]
  res <- evaluate_pipeline(fit_effect(null_expr))
  list(iter = i, n = res$selected, r = res$r,
       genes = res$genes, a = res$a_lfc, b = res$b_lfc,
       elig = res$eligible, bsig = res$b_significant)
}

t0 <- Sys.time()
if (nworkers > 1L) {
  cl <- makeCluster(nworkers)
  clusterEvalQ(cl, suppressPackageStartupMessages(library(limma)))
  clusterExport(cl, c("permutations", "nuisance_fitted", "nuisance_residuals",
                      "design_full", "fit_effect", "evaluate_pipeline",
                      "anchor", "hubs", "y_chromosome_genes",
                      "fdr_cut", "lfc_cut", "min_selected"),
                envir = environment())
  results <- parLapply(cl, idx, worker)
  stopCluster(cl)
} else {
  results <- lapply(idx, worker)
}
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("chunk %d-%d: %.2f min (%.2f s/draw)\n",
            chunk_start, max(idx), elapsed, 60 * elapsed / length(idx)))

gene_rows <- do.call(rbind, lapply(results, function(z)
  data.frame(Iteration = z$iter, Gene = z$genes, A_log2FC = z$a, B_log2FC = z$b,
             stringsAsFactors = FALSE)))
iter_rows <- do.call(rbind, lapply(results, function(z)
  data.frame(Iteration = z$iter, SignificantTestArmGenes = z$bsig,
             EligiblePoolGenes = z$elig, SameSignSelectedGenes = z$n,
             SelectedPearson = z$r, stringsAsFactors = FALSE)))

write_csv_atomic(gene_rows, p("outputs", "positive_control_masld_masld",
  sprintf("replay_genes_%04d_%04d.csv", chunk_start, max(idx))))
write_csv_atomic(iter_rows, p("outputs", "positive_control_masld_masld",
  sprintf("replay_iterations_%04d_%04d.csv", chunk_start, max(idx))))
cat("chunk saved.\n")

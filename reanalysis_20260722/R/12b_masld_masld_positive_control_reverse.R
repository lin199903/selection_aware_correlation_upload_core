script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE)), "_common.R"))

required_packages(c("limma", "jsonlite"))
suppressPackageStartupMessages(library(limma))

## 12b — Reverse-direction sensitivity for the same-disease MASLD stress test:
## fixed arm = GSE135251 (limma vst effects), replay arm = GSE126848 (limma on vst, ~group).
## Usage: Rscript 12b_..._reverse.R obs | chunk <start> <size>

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
stopifnot(mode %in% c("obs", "chunk"))

spec <- read_spec()
ensure_dirs(p("outputs", "positive_control_masld_masld_reverse"), p("logs"))

seed <- as.integer(spec$seed_benchmark) + 32000L   # 20290725, frozen for reverse
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

## fixed arm: GSE135251 limma effects (identical estimator to primary test arm)
g251 <- readRDS(normalizePath(p("data", "derived", "GSE135251_rebuilt.rds"), winslash = "/", mustWork = TRUE))
expr251 <- g251$vst
meta251 <- g251$metadata
meta251$group <- factor(meta251$group, levels = c("Control", "MASLD"))
design251 <- model.matrix(~ group, data = meta251)
fit251 <- eBayes(lmFit(expr251, design251), trend = TRUE, robust = TRUE)
tt251 <- topTable(fit251, coef = "groupMASLD", number = Inf, sort.by = "none")
anchor <- data.frame(Gene = rownames(tt251), A_log2FC = as.numeric(tt251$logFC),
                     A_FDR = as.numeric(tt251$adj.P.Val), stringsAsFactors = FALSE)
anchor$A_significant <- is.finite(anchor$A_FDR) & anchor$A_FDR < fdr_cut &
  is.finite(anchor$A_log2FC) & abs(anchor$A_log2FC) > lfc_cut

hubs <- unique(trimws(read.csv(
  p("outputs", "network", "MASLD_topology_hubs_primary.csv"), check.names = FALSE)$Gene))
hubs <- hubs[nzchar(hubs)]

## replay arm: GSE126848
g126 <- readRDS(p("data", "derived", "GSE126848_rebuilt.rds"))
expr <- g126$vst
coldata <- g126$coldata
cat("coldata cols:", paste(names(coldata), collapse=", "), "\n")
grp_col <- names(coldata)[1]
coldata$group <- factor(coldata[[grp_col]])
cat("group levels:", paste(levels(coldata$group), collapse=", "), "\n")
ctrl_level <- levels(coldata$group)[1]
design_full <- model.matrix(~ group, data = coldata)
design_nuisance <- model.matrix(~ 1, data = coldata)
fit_nuisance <- lmFit(expr, design_nuisance)
nuisance_fitted <- fit_nuisance$coefficients %*% t(design_nuisance)
nuisance_residuals <- expr - nuisance_fitted
coef_name <- colnames(design_full)[2]
cat("effect coefficient:", coef_name, "\n")

fit_effect <- function(m) {
  fit <- eBayes(lmFit(m, design_full), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
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
  write_csv_atomic(data.frame(Gene = obs$genes, A_log2FC = obs$a_lfc, B_log2FC = obs$b_lfc),
                   p("outputs", "positive_control_masld_masld_reverse", "observed_selected_set.csv"))
  write_csv_atomic(data.frame(
    Analysis = "MASLD_MASLD_positive_control_reverse_observed",
    CommonGenes = obs$common, TestArmSignificant = obs$b_significant,
    EligiblePoolGenes = obs$eligible, SameSignSelectedGenes = obs$selected,
    ObservedSelectedPearson = obs$r),
    p("outputs", "positive_control_masld_masld_reverse", "observed_summary.csv"))
  quit(save = "no")
}

chunk_start <- as.integer(args[2]); chunk_size <- as.integer(args[3])
set.seed(seed)
permutations <- replicate(n_perm, sample.int(ncol(expr)), simplify = FALSE)
idx <- chunk_start:min(chunk_start + chunk_size - 1L, n_perm)
rows <- vector("list", length(idx))
t0 <- Sys.time()
for (k in seq_along(idx)) {
  i <- idx[k]
  null_expr <- nuisance_fitted + nuisance_residuals[, permutations[[i]], drop = FALSE]
  res <- evaluate_pipeline(fit_effect(null_expr))
  rows[[k]] <- data.frame(Iteration = i, SignificantTestArmGenes = res$b_significant,
                          EligiblePoolGenes = res$eligible, SameSignSelectedGenes = res$selected,
                          SelectedPearson = res$r, stringsAsFactors = FALSE)
}
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("chunk %d-%d: %.2f min (%.2f s/draw)\n", chunk_start, max(idx), elapsed,
            60 * elapsed / length(idx)))
write_csv_atomic(do.call(rbind, rows), p("outputs", "positive_control_masld_masld_reverse",
  sprintf("replay_iterations_%04d_%04d.csv", chunk_start, max(idx))))
cat("chunk saved.\n")

script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE)), "_common.R"))

required_packages(c("limma", "matrixStats", "jsonlite"))
suppressPackageStartupMessages({
  library(limma)
  library(matrixStats)
})

spec <- read_spec()
set.seed(spec$seed_primary + 10L)
ensure_dirs(p("outputs", "osa_cpap"), p("logs"))

input_path <- p("data", "derived", "GSE135917_rebuilt.rds")
if (!file.exists(input_path)) stop("Run 08_osa_de.R first")
x <- readRDS(input_path)
expr <- x$gene_expression_all_samples
metadata <- x$cpap_metadata
metadata <- metadata[match(colnames(expr), metadata$sample_id, nomatch = 0L), , drop = FALSE]
metadata <- metadata[metadata$study_group == "STUDY GROUP 2", , drop = FALSE]
expr <- expr[, metadata$sample_id, drop = FALSE]
stopifnot(ncol(expr) == 48L, !anyDuplicated(metadata$sample_id))

metadata$condition <- factor(
  ifelse(metadata$group == "CPAP_post", "Post", "Pre"),
  levels = c("Pre", "Post")
)
metadata$patient_id <- factor(metadata$patient_id)
pair_table <- table(metadata$patient_id, metadata$condition)
stopifnot(nrow(pair_table) == 24L, all(pair_table == 1L))

# Paired limma model. Patient fixed effects absorb all stable person-level
# characteristics; conditionPost estimates the within-patient CPAP change.
design <- model.matrix(~ patient_id + condition, data = metadata)
stopifnot(qr(design)$rank == ncol(design), "conditionPost" %in% colnames(design))
fit <- eBayes(lmFit(expr, design), trend = TRUE, robust = TRUE)
de <- topTable(fit, coef = "conditionPost", number = Inf, sort.by = "none")
de$Gene <- rownames(de)
de$Significant <- is.finite(de$adj.P.Val) &
  de$adj.P.Val < as.numeric(spec$de_thresholds$fdr) &
  abs(de$logFC) > as.numeric(spec$de_thresholds$abs_log2fc)
de <- de[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "Significant")]
write_csv_atomic(de, p("outputs", "osa_cpap", "GSE135917_CPAP_paired_limma.csv"))

# These compact axes are fixed biological summaries, not selected from the
# CPAP results. Scores are the mean of gene-wise z scores across the 48 paired
# samples. The remodelling coordinate follows its signed definition.
gene_sets <- list(
  IFN_STAT1 = c("STAT1", "IRF1", "IRF7", "IFIT1", "IFIT2", "IFIT3", "IFI44",
                "IFITM3", "MX1", "OAS2", "OAS3", "OASL", "DDX60", "CXCL9", "CXCL10", "CXCL11"),
  Inflammasome_TLR = c("NLRP3", "TLR2", "TLR4", "IL1B", "CASP1", "PYCARD",
                       "PTGS2", "CCL2", "CCL20", "CXCL2", "S100A8", "S100A9"),
  TNF_NFKB_feedback = c("NFKBIA", "TNFAIP3", "RELB", "BCL3", "SOCS1", "SOCS3",
                        "TNFAIP2", "ICAM1", "CXCL8"),
  AP1_core = c("FOS", "JUN", "JUNB", "JUND")
)

z_rows <- function(m) {
  out <- t(apply(m, 1, function(v) as.numeric(scale(v))))
  rownames(out) <- rownames(m)
  colnames(out) <- colnames(m)
  out
}
scores <- list()
coverage <- list()
for (nm in names(gene_sets)) {
  genes <- intersect(gene_sets[[nm]], rownames(expr))
  coverage[[nm]] <- data.frame(
    Axis = nm, RequestedGenes = length(gene_sets[[nm]]),
    AvailableGenes = length(genes), Genes = paste(genes, collapse = ";"),
    stringsAsFactors = FALSE
  )
  scores[[nm]] <- colMeans(z_rows(expr[genes, , drop = FALSE]), na.rm = TRUE)
}
remodel_genes <- c("ATF3", "JUNB", "MAFB")
stopifnot(all(remodel_genes %in% rownames(expr)))
remodel_z <- z_rows(expr[remodel_genes, , drop = FALSE])
scores$Remodel_ATF3_high_JUNB_MAFB_low <- remodel_z["ATF3", ] -
  remodel_z["JUNB", ] - remodel_z["MAFB", ]
coverage[["Remodel_ATF3_high_JUNB_MAFB_low"]] <- data.frame(
  Axis = "Remodel_ATF3_high_JUNB_MAFB_low", RequestedGenes = 3L,
  AvailableGenes = 3L, Genes = paste(remodel_genes, collapse = ";"),
  stringsAsFactors = FALSE
)

score_matrix <- do.call(cbind, scores)
score_long <- do.call(rbind, lapply(seq_len(ncol(score_matrix)), function(j) {
  data.frame(
    sample_id = rownames(score_matrix), patient_id = as.character(metadata$patient_id),
    condition = as.character(metadata$condition), Axis = colnames(score_matrix)[j],
    Score = score_matrix[, j], stringsAsFactors = FALSE
  )
}))

axis_results <- do.call(rbind, lapply(colnames(score_matrix), function(axis) {
  pre_idx <- metadata$condition == "Pre"
  post_idx <- metadata$condition == "Post"
  pre <- score_matrix[pre_idx, axis]
  post <- score_matrix[post_idx, axis]
  pre_id <- as.character(metadata$patient_id[pre_idx])
  post_id <- as.character(metadata$patient_id[post_idx])
  post <- post[match(pre_id, post_id)]
  delta <- post - pre
  tt <- t.test(delta, mu = 0)
  wt <- suppressWarnings(wilcox.test(delta, mu = 0, exact = FALSE))
  data.frame(
    Axis = axis, Pairs = length(delta), MeanPre = mean(pre), MeanPost = mean(post),
    MeanDelta = mean(delta), SEDelta = sd(delta) / sqrt(length(delta)),
    T = unname(tt$statistic), PairedT_P = tt$p.value,
    WilcoxonV = unname(wt$statistic), PairedWilcoxon_P = wt$p.value,
    stringsAsFactors = FALSE
  )
}))
axis_results$PairedT_FDR <- p.adjust(axis_results$PairedT_P, "BH")
axis_results$PairedWilcoxon_FDR <- p.adjust(axis_results$PairedWilcoxon_P, "BH")

write_csv_atomic(do.call(rbind, coverage), p("outputs", "osa_cpap", "CPAP_axis_gene_coverage.csv"))
write_csv_atomic(score_long, p("outputs", "osa_cpap", "CPAP_axis_scores_long.csv"))
write_csv_atomic(axis_results, p("outputs", "osa_cpap", "CPAP_axis_paired_results.csv"))

gene_check <- de[match(c("NFKBIA", "RELB", "STAT1", "IRF1", "JUNB", "ATF3", "MAFB"), de$Gene), , drop = FALSE]
write_csv_atomic(gene_check, p("outputs", "osa_cpap", "CPAP_prespecified_gene_results.csv"))

plot_axes <- function(device, file) {
  if (device == "pdf") pdf(file, width = 8, height = 5.5)
  if (device == "png") png(file, width = 2400, height = 1650, res = 300)
  on.exit(dev.off(), add = TRUE)
  ord <- order(axis_results$MeanDelta)
  y <- seq_along(ord)
  lo <- axis_results$MeanDelta[ord] - 1.96 * axis_results$SEDelta[ord]
  hi <- axis_results$MeanDelta[ord] + 1.96 * axis_results$SEDelta[ord]
  plot(axis_results$MeanDelta[ord], y, xlim = range(c(lo, hi, 0)), pch = 19,
       yaxt = "n", ylab = "", xlab = "Post-CPAP minus pre-CPAP score (95% CI)",
       main = "GSE135917 paired intervention")
  segments(lo, y, hi, y, lwd = 2, col = "#2166AC")
  points(axis_results$MeanDelta[ord], y, pch = 19, col = "#B2182B")
  axis(2, at = y, labels = axis_results$Axis[ord], las = 1, cex.axis = 0.8)
  abline(v = 0, lty = 2, col = "grey40")
}
plot_axes("pdf", p("outputs", "osa_cpap", "CPAP_axis_changes.pdf"))
plot_axes("png", p("outputs", "osa_cpap", "CPAP_axis_changes.png"))

saveRDS(list(
  spec_version = spec$spec_version, metadata = metadata, design = design,
  paired_de = de, gene_sets = gene_sets, score_matrix = score_matrix,
  axis_results = axis_results
), p("outputs", "osa_cpap", "GSE135917_CPAP_rebuilt.rds"), compress = "xz")

contract <- list(
  created = format(Sys.time(), tz = "UTC", usetz = TRUE),
  spec_version = spec$spec_version,
  accession = "GSE135917",
  population = "Study Group 2 paired CPAP cohort",
  n_pairs = nrow(pair_table),
  design_rank = qr(design)$rank,
  design_columns = ncol(design),
  n_gene_level_significant = sum(de$Significant),
  n_axes = nrow(axis_results),
  n_axis_t_fdr_lt_0_05 = sum(axis_results$PairedT_FDR < 0.05),
  n_axis_wilcoxon_fdr_lt_0_05 = sum(axis_results$PairedWilcoxon_FDR < 0.05),
  paired_design = TRUE
)
stopifnot(contract$n_pairs == 24L, contract$design_rank == 25L,
          contract$n_axes == 5L)
write_json_atomic(contract, p("outputs", "osa_cpap", "osa_cpap_contract.json"))
save_session_info("10_osa_cpap_intervention")
print(axis_results)
print(gene_check)
cat("OSA paired CPAP intervention rebuild completed.\n")

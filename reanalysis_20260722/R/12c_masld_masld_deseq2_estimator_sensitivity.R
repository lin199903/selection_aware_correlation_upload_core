script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE)), "_common.R"))

required_packages(c("DESeq2", "limma", "jsonlite"))
suppressPackageStartupMessages({library(DESeq2); library(limma)})

## Observed-effect estimator sensitivity for the GSE135251 same-disease stress test.
## This does NOT replace the 2,000-draw limma replay. It asks whether replacing
## the observed GSE135251 effect estimator by DESeq2 materially changes the
## selected set or selected correlation under the same eligibility rule.
##
## Required author-side object: data/derived/GSE135251_rebuilt.rds
## Expected components: metadata, vst, and a raw integer count matrix under one
## of: counts, raw_counts, count_matrix. Raw GEO HTSeq counts are third-party
## inputs and are not redistributed in the public release.

spec <- read_spec()
ensure_dirs(p("outputs", "positive_control_masld_masld", "estimator_sensitivity"), p("logs"))
fdr_cut <- as.numeric(spec$de_thresholds$fdr)
lfc_cut <- as.numeric(spec$de_thresholds$abs_log2fc)
min_selected <- as.integer(spec$rule_replay$min_selected_genes)

y_chromosome_genes <- c(
  "UTY", "USP9Y", "DDX3Y", "KDM5D", "EIF1AY", "ZFY", "SRY",
  "NLGN4Y", "RPS4Y1", "RPS4Y2", "TSPY1", "RBMY1A1", "DAZ1",
  "PRKY", "AMELY", "TBL1Y", "PCDH11Y", "TMSB4Y", "VCY",
  "CDY1", "CDY2A", "HSFY1", "TXLNGY", "BPY2", "PRY"
)

obj_path <- p("data", "derived", "GSE135251_rebuilt.rds")
if (!file.exists(obj_path)) {
  stop("Missing author-side GSE135251_rebuilt.rds. Reconstruct the public GEO HTSeq counts before running this sensitivity.")
}
g251 <- readRDS(obj_path)
if (is.null(g251$metadata) || is.null(g251$vst)) stop("GSE135251 object lacks metadata/vst components")
count_candidates <- c("counts", "raw_counts", "count_matrix")
count_name <- count_candidates[vapply(count_candidates, function(z) !is.null(g251[[z]]), logical(1))]
if (!length(count_name)) {
  stop("GSE135251 object has no raw-count component (counts/raw_counts/count_matrix). DESeq2 sensitivity requires the public HTSeq count matrix.")
}
counts <- as.matrix(g251[[count_name[1]]])
storage.mode(counts) <- "integer"
meta <- g251$metadata
if (!"sample_id" %in% names(meta)) stop("metadata lacks sample_id")
if (!identical(colnames(counts), meta$sample_id)) counts <- counts[, meta$sample_id, drop = FALSE]
if (!"group" %in% names(meta)) stop("metadata lacks group")
meta$group <- factor(meta$group, levels = c("Control", "MASLD"))
if (anyNA(meta$group)) stop("group must map to Control/MASLD")

## Match the light count filter used elsewhere in the MASLD reconstruction.
keep <- rowSums(counts >= 10L) >= 6L
counts <- counts[keep, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ group)
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("group", "MASLD", "Control"), independentFiltering = TRUE)
deseq <- data.frame(
  Gene = rownames(res),
  B_log2FC_DESeq2 = as.numeric(res$log2FoldChange),
  B_FDR_DESeq2 = as.numeric(res$padj),
  stringsAsFactors = FALSE
)
deseq <- deseq[nzchar(deseq$Gene) & !duplicated(deseq$Gene), , drop = FALSE]
deseq$B_significant_DESeq2 <- is.finite(deseq$B_FDR_DESeq2) & deseq$B_FDR_DESeq2 < fdr_cut &
  is.finite(deseq$B_log2FC_DESeq2) & abs(deseq$B_log2FC_DESeq2) > lfc_cut

## Refit the frozen limma observed estimator from the same author-side object.
expr <- g251$vst
if (!identical(colnames(expr), meta$sample_id)) expr <- expr[, meta$sample_id, drop = FALSE]
design <- model.matrix(~ group, data = meta)
fit <- eBayes(lmFit(expr, design), trend = TRUE, robust = TRUE)
tt <- topTable(fit, coef = "groupMASLD", number = Inf, sort.by = "none")
limma_eff <- data.frame(
  Gene = rownames(tt),
  B_log2FC_limma = as.numeric(tt$logFC),
  B_FDR_limma = as.numeric(tt$adj.P.Val),
  stringsAsFactors = FALSE
)
limma_eff <- limma_eff[nzchar(limma_eff$Gene) & !duplicated(limma_eff$Gene), , drop = FALSE]

anchor_de <- read.csv(p("outputs", "de", "GSE126848_DESeq2_primary.csv"), check.names = FALSE)
anchor <- data.frame(Gene = trimws(as.character(anchor_de$Gene)),
                     A_log2FC = as.numeric(anchor_de$log2FoldChange),
                     A_FDR = as.numeric(anchor_de$padj), stringsAsFactors = FALSE)
anchor <- anchor[nzchar(anchor$Gene) & !duplicated(anchor$Gene), , drop = FALSE]
anchor$A_significant <- is.finite(anchor$A_FDR) & anchor$A_FDR < fdr_cut &
  is.finite(anchor$A_log2FC) & abs(anchor$A_log2FC) > lfc_cut
hubs <- unique(trimws(read.csv(p("outputs", "network", "MASLD_topology_hubs_primary.csv"), check.names = FALSE)$Gene))
hubs <- hubs[nzchar(hubs)]

calc_selected <- function(b, lfc_col, fdr_col, label) {
  b$B_log2FC <- b[[lfc_col]]; b$B_FDR <- b[[fdr_col]]
  b$B_significant <- is.finite(b$B_FDR) & b$B_FDR < fdr_cut & is.finite(b$B_log2FC) & abs(b$B_log2FC) > lfc_cut
  d <- merge(anchor, b[, c("Gene", "B_log2FC", "B_FDR", "B_significant")], by = "Gene", all = FALSE, sort = FALSE)
  d <- d[is.finite(d$A_log2FC) & is.finite(d$B_log2FC) & !d$Gene %in% y_chromosome_genes, , drop = FALSE]
  d$eligible <- d$Gene %in% hubs & (d$A_significant | d$B_significant)
  pool <- d[d$eligible, , drop = FALSE]
  sel <- pool[pool$A_log2FC * pool$B_log2FC > 0, , drop = FALSE]
  r <- if (nrow(sel) >= min_selected) cor(sel$A_log2FC, sel$B_log2FC) else NA_real_
  list(label = label, all = d, pool = pool, selected = sel, r = r)
}

s_deseq <- calc_selected(deseq, "B_log2FC_DESeq2", "B_FDR_DESeq2", "DESeq2")
s_limma <- calc_selected(limma_eff, "B_log2FC_limma", "B_FDR_limma", "limma")

cmp <- merge(deseq, limma_eff, by = "Gene", all = FALSE)
cmp <- cmp[is.finite(cmp$B_log2FC_DESeq2) & is.finite(cmp$B_log2FC_limma), , drop = FALSE]
pearson <- cor(cmp$B_log2FC_DESeq2, cmp$B_log2FC_limma)
spearman <- cor(cmp$B_log2FC_DESeq2, cmp$B_log2FC_limma, method = "spearman")
sign_agree <- mean(sign(cmp$B_log2FC_DESeq2) == sign(cmp$B_log2FC_limma))
inter <- length(intersect(s_deseq$selected$Gene, s_limma$selected$Gene))
uni <- length(union(s_deseq$selected$Gene, s_limma$selected$Gene))
jaccard <- if (uni) inter / uni else NA_real_

summary <- data.frame(
  Analysis = "GSE135251_DESeq2_vs_limma_observed_estimator_sensitivity",
  CommonEffectGenes = nrow(cmp),
  EffectPearson = pearson,
  EffectSpearman = spearman,
  SignAgreement = sign_agree,
  LimmaEligible = nrow(s_limma$pool),
  DESeq2Eligible = nrow(s_deseq$pool),
  LimmaSelected = nrow(s_limma$selected),
  DESeq2Selected = nrow(s_deseq$selected),
  SelectedSetJaccard = jaccard,
  LimmaSelectedPearson = s_limma$r,
  DESeq2SelectedPearson = s_deseq$r,
  stringsAsFactors = FALSE
)
write_csv_atomic(summary, p("outputs", "positive_control_masld_masld", "estimator_sensitivity", "deseq2_vs_limma_summary.csv"))
write_csv_atomic(cmp, p("outputs", "positive_control_masld_masld", "estimator_sensitivity", "deseq2_vs_limma_gene_effects.csv"))
write_csv_atomic(s_deseq$selected, p("outputs", "positive_control_masld_masld", "estimator_sensitivity", "deseq2_selected_set.csv"))
print(summary)
save_session_info("12c_masld_masld_deseq2_estimator_sensitivity")

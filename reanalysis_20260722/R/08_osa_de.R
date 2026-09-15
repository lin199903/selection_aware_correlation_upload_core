script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE)), "_common.R"))

required_packages(c(
  "GEOquery", "Biobase", "limma", "AnnotationDbi",
  "hugene10sttranscriptcluster.db", "matrixStats", "jsonlite"
))
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(limma)
  library(AnnotationDbi)
  library(hugene10sttranscriptcluster.db)
  library(matrixStats)
})

spec <- read_spec()
set.seed(spec$seed_primary)
ensure_dirs(p("outputs", "osa_de"), p("logs"))

series_path <- p("data", "raw", "GSE135917_series_matrix.txt.gz")
metadata_path <- p("data", "derived", "metadata_GSE135917.csv")
if (!file.exists(series_path) || !file.exists(metadata_path)) stop("Run 00_audit_inputs.R first")

eset <- getGEO(filename = series_path, getGPL = FALSE)
if (is.list(eset)) eset <- eset[[1L]]
probe_expr <- exprs(eset)
metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
metadata <- metadata[match(colnames(probe_expr), metadata$sample_id), , drop = FALSE]
stopifnot(!anyNA(metadata$sample_id), identical(colnames(probe_expr), metadata$sample_id))
stopifnot(all(metadata$platform == "GPL6244"), nrow(metadata) == 66L)

discovery <- metadata$study_group == "STUDY GROUP 1"
stopifnot(sum(discovery) == 18L,
          sum(metadata$group[discovery] == "Control") == 8L,
          sum(metadata$group[discovery] == "OSA_baseline") == 10L)

# GPL6244 is the Affymetrix Human Gene 1.0 ST transcript-cluster platform.
# Probes mapping to more than one current HGNC symbol are removed. When several
# unambiguous probes map to the same symbol, the probe with the highest MAD in
# the fixed Study Group 1 discovery cohort is retained.
annotation <- AnnotationDbi::select(
  hugene10sttranscriptcluster.db,
  keys = rownames(probe_expr), keytype = "PROBEID",
  columns = c("SYMBOL", "ENTREZID", "GENENAME")
)
annotation <- annotation[!is.na(annotation$SYMBOL) & nzchar(annotation$SYMBOL), , drop = FALSE]
probe_multiplicity <- table(annotation$PROBEID)
ambiguous_probes <- names(probe_multiplicity)[probe_multiplicity > 1L]
annotation <- annotation[!annotation$PROBEID %in% ambiguous_probes, , drop = FALSE]
annotation <- annotation[!duplicated(annotation$PROBEID), , drop = FALSE]
annotation$DiscoveryMAD <- rowMads(
  probe_expr[match(annotation$PROBEID, rownames(probe_expr)), discovery, drop = FALSE]
)
annotation <- annotation[is.finite(annotation$DiscoveryMAD), , drop = FALSE]
annotation <- annotation[order(annotation$SYMBOL, -annotation$DiscoveryMAD, annotation$PROBEID), , drop = FALSE]
selected_annotation <- annotation[!duplicated(annotation$SYMBOL), , drop = FALSE]

gene_expr <- probe_expr[match(selected_annotation$PROBEID, rownames(probe_expr)), , drop = FALSE]
rownames(gene_expr) <- selected_annotation$SYMBOL
storage.mode(gene_expr) <- "numeric"
stopifnot(!anyDuplicated(rownames(gene_expr)), !anyNA(gene_expr))

coldata <- metadata[discovery, , drop = FALSE]
rownames(coldata) <- coldata$sample_id
expr_discovery <- gene_expr[, coldata$sample_id, drop = FALSE]
coldata$group <- factor(
  ifelse(coldata$group == "Control", "Control", "OSA"),
  levels = c("Control", "OSA")
)
coldata$sex <- factor(coldata$sex)
coldata$age_z <- as.numeric(scale(coldata$age))
coldata$bmi_z <- as.numeric(scale(coldata$bmi))

design_primary <- model.matrix(~ age_z + bmi_z + sex + group, data = coldata)
design_sensitivity <- model.matrix(~ group, data = coldata)
stopifnot(qr(design_primary)$rank == ncol(design_primary),
          "groupOSA" %in% colnames(design_primary))

fit_limma <- function(design) {
  eBayes(lmFit(expr_discovery, design), trend = TRUE, robust = TRUE)
}
format_de <- function(fit, design_label) {
  out <- topTable(fit, coef = "groupOSA", number = Inf, sort.by = "none")
  out$Gene <- rownames(out)
  out$ProbeID <- selected_annotation$PROBEID[match(out$Gene, selected_annotation$SYMBOL)]
  out$Design <- design_label
  out$Significant <- is.finite(out$adj.P.Val) &
    out$adj.P.Val < as.numeric(spec$de_thresholds$fdr) &
    abs(out$logFC) > as.numeric(spec$de_thresholds$abs_log2fc)
  out[, c("Gene", "ProbeID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "Design", "Significant")]
}

de_primary <- format_de(fit_limma(design_primary), "age_z + bmi_z + sex + group")
de_sensitivity <- format_de(fit_limma(design_sensitivity), "group")

write_csv_atomic(selected_annotation, p("outputs", "osa_de", "GSE135917_probe_to_gene_primary.csv"))
write_csv_atomic(coldata, p("outputs", "osa_de", "GSE135917_discovery_metadata.csv"), row.names = FALSE)
write_csv_atomic(de_primary, p("outputs", "osa_de", "GSE135917_limma_primary.csv"))
write_csv_atomic(de_sensitivity, p("outputs", "osa_de", "GSE135917_limma_unadjusted_sensitivity.csv"))

cpap <- metadata[metadata$study_group == "STUDY GROUP 2", , drop = FALSE]
pair_table <- table(cpap$patient_id, cpap$group)
stopifnot(nrow(pair_table) == 24L,
          all(pair_table[, "OSA_baseline"] == 1L),
          all(pair_table[, "CPAP_post"] == 1L))
write_csv_atomic(cpap, p("outputs", "osa_de", "GSE135917_CPAP_metadata.csv"))

plot_pca <- function(device, file) {
  if (device == "pdf") pdf(file, width = 7, height = 5.5)
  if (device == "png") png(file, width = 2100, height = 1650, res = 300)
  on.exit(dev.off(), add = TRUE)
  mad <- rowMads(expr_discovery)
  top <- order(mad, decreasing = TRUE)[seq_len(min(1000L, length(mad)))]
  pc <- prcomp(t(expr_discovery[top, , drop = FALSE]), scale. = TRUE)
  cols <- c(Control = "#2166AC", OSA = "#B2182B")
  plot(pc$x[, 1], pc$x[, 2], pch = 19, cex = 1.3,
       col = cols[as.character(coldata$group)],
       xlab = sprintf("PC1 (%.1f%%)", 100 * summary(pc)$importance[2, 1]),
       ylab = sprintf("PC2 (%.1f%%)", 100 * summary(pc)$importance[2, 2]),
       main = "GSE135917 Study Group 1")
  legend("topright", legend = names(cols), col = cols, pch = 19, bty = "n")
}
plot_pca("pdf", p("outputs", "osa_de", "GSE135917_discovery_PCA.pdf"))
plot_pca("png", p("outputs", "osa_de", "GSE135917_discovery_PCA.png"))

saveRDS(list(
  spec_version = spec$spec_version,
  platform = "GPL6244",
  processing = unique(pData(eset)$data_processing),
  gene_expression_all_samples = gene_expr,
  discovery_expression = expr_discovery,
  discovery_coldata = coldata,
  cpap_metadata = cpap,
  selected_annotation = selected_annotation,
  de_primary = de_primary,
  de_unadjusted_sensitivity = de_sensitivity,
  design_primary = design_primary,
  design_sensitivity = design_sensitivity
), p("data", "derived", "GSE135917_rebuilt.rds"), compress = "xz")

contract <- list(
  created = format(Sys.time(), tz = "UTC", usetz = TRUE),
  spec_version = spec$spec_version,
  accession = "GSE135917",
  platform = "GPL6244",
  source_processing = unique(pData(eset)$data_processing),
  n_series_samples = ncol(probe_expr),
  n_discovery_samples = ncol(expr_discovery),
  discovery_groups = as.list(table(coldata$group)),
  discovery_sex = as.list(table(coldata$sex)),
  n_cpap_pairs = nrow(pair_table),
  n_platform_probe_rows = nrow(probe_expr),
  n_ambiguous_probes_removed = length(ambiguous_probes),
  n_gene_symbols = nrow(gene_expr),
  primary_design_rank = qr(design_primary)$rank,
  primary_design_columns = colnames(design_primary),
  n_primary_significant = sum(de_primary$Significant),
  n_unadjusted_significant = sum(de_sensitivity$Significant),
  primary_min_fdr = min(de_primary$adj.P.Val, na.rm = TRUE),
  probe_collapse_rule = "drop multi-symbol probes; highest discovery-cohort MAD probe per symbol"
)
stopifnot(contract$platform == "GPL6244",
          contract$n_discovery_samples == 18L,
          contract$n_cpap_pairs == 24L,
          contract$n_gene_symbols > 15000L,
          contract$primary_design_rank == 5L)
write_json_atomic(contract, p("outputs", "osa_de", "osa_de_contract.json"))
save_session_info("08_osa_de")
cat("OSA source-level differential-expression rebuild completed.\n")

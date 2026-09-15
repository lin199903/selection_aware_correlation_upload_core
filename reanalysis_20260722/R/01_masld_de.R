script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE)), "_common.R"))

required_packages(c(
  "DESeq2", "GEOquery", "Biobase", "data.table", "AnnotationDbi",
  "org.Hs.eg.db", "jsonlite", "ggplot2"
))

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggplot2)
})

spec <- read_spec()
set.seed(spec$seed_primary)
ensure_dirs(p("data", "derived"), p("outputs", "de"), p("logs"))

contract_path <- p("data", "derived", "data_contract.json")
if (!file.exists(contract_path)) stop("Run 00_audit_inputs.R first")

meta126 <- read.csv(p("data", "derived", "metadata_GSE126848.csv"), check.names = FALSE)
meta130 <- read.csv(p("data", "derived", "metadata_GSE130970.csv"), check.names = FALSE)

read_126848_counts <- function() {
  x <- fread(p("data", "raw", "GSE126848_Gene_counts_raw.txt.gz"),
             header = FALSE, showProgress = FALSE)
  header <- as.character(unlist(x[1L, -1L, with = FALSE], use.names = FALSE))
  gene_ids <- sub("\\..*$", "", as.character(x[-1L, 1L][[1L]]))
  mat <- as.matrix(x[-1L, -1L, with = FALSE])
  storage.mode(mat) <- "numeric"
  stopifnot(ncol(mat) == nrow(meta126), identical(header, as.character(meta126$raw_column_id)))
  colnames(mat) <- meta126$sample_id
  symbols <- mapIds(org.Hs.eg.db, keys = gene_ids, column = "SYMBOL",
                    keytype = "ENSEMBL", multiVals = "first")
  collapse_counts_by_symbol(mat, symbols)
}

read_130970_counts <- function() {
  x <- fread(p("data", "raw", "GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv.gz"),
             showProgress = FALSE)
  ids <- as.character(x[[1L]])
  mat <- as.matrix(x[, -1L, with = FALSE])
  storage.mode(mat) <- "numeric"
  stopifnot(identical(colnames(mat), meta130$title))
  colnames(mat) <- meta130$sample_id
  symbols <- mapIds(org.Hs.eg.db, keys = ids, column = "SYMBOL",
                    keytype = "ENTREZID", multiVals = "first")
  collapse_counts_by_symbol(mat, symbols)
}

counts126_all <- read_126848_counts()
counts130 <- read_130970_counts()
stopifnot(!anyDuplicated(rownames(counts126_all)), !anyDuplicated(rownames(counts130)))

filter_counts <- function(counts, min_count, min_samples) {
  counts[rowSums(counts >= min_count) >= min_samples, , drop = FALSE]
}
counts126_all <- filter_counts(counts126_all, spec$masld_primary$min_count,
                               spec$masld_primary$min_samples)
counts130 <- filter_counts(counts130, spec$masld_replication$min_count,
                           spec$masld_replication$min_samples)

col126_all <- data.frame(
  group = factor(meta126$group, levels = c("Control", "MASLD")),
  sex = factor(toupper(substr(meta126$sex, 1L, 1L))),
  row.names = meta126$sample_id
)
male_ids <- rownames(col126_all)[col126_all$sex == "M"]
col126 <- col126_all[male_ids, "group", drop = FALSE]
counts126 <- counts126_all[, male_ids, drop = FALSE]
counts126 <- filter_counts(counts126, spec$masld_primary$min_count,
                           spec$masld_primary$min_samples)
col130 <- data.frame(
  group = factor(meta130$group, levels = c("Control", "MASLD")),
  sex = factor(toupper(substr(meta130$sex, 1L, 1L))),
  age_z = as.numeric(scale(meta130$age)),
  row.names = meta130$sample_id
)
stopifnot(identical(colnames(counts126), rownames(col126)))
stopifnot(identical(colnames(counts126_all), rownames(col126_all)))
stopifnot(identical(colnames(counts130), rownames(col130)))

fit_deseq <- function(counts, coldata, design, contrast = c("group", "MASLD", "Control")) {
  rounded <- round(counts)
  storage.mode(rounded) <- "integer"
  mm <- model.matrix(design, coldata)
  if (qr(mm)$rank != ncol(mm)) stop("Rank-deficient design: ", deparse(design))
  dds <- DESeqDataSetFromMatrix(rounded, coldata, design = design)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = contrast, alpha = spec$de_thresholds$fdr)
  vst_obj <- vst(dds, blind = FALSE)
  list(dds = dds, result = res, vst = assay(vst_obj), design_rank = qr(mm)$rank,
       design_columns = colnames(mm))
}

fit126 <- fit_deseq(counts126, col126, ~ group)
fit126_full <- fit_deseq(counts126_all, col126_all, ~ sex + group)
fit130 <- fit_deseq(counts130, col130, ~ sex + age_z + group)
fit130_unadjusted <- fit_deseq(counts130, col130[, "group", drop = FALSE], ~ group)

result_table <- function(res) {
  z <- as.data.frame(res)
  z$Gene <- rownames(z)
  z <- z[, c("Gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  z$significant <- !is.na(z$padj) & z$padj < spec$de_thresholds$fdr &
    abs(z$log2FoldChange) > spec$de_thresholds$abs_log2fc
  z[order(z$padj, -abs(z$log2FoldChange), na.last = TRUE), ]
}

de126 <- result_table(fit126$result)
de126_full <- result_table(fit126_full$result)
de130 <- result_table(fit130$result)
de130_unadj <- result_table(fit130_unadjusted$result)
write_csv_atomic(de126, p("outputs", "de", "GSE126848_DESeq2_primary.csv"))
write_csv_atomic(de126[de126$significant, ], p("outputs", "de", "GSE126848_DESeq2_significant.csv"))
write_csv_atomic(de126_full, p("outputs", "de", "GSE126848_DESeq2_full_sex_adjusted_sensitivity.csv"))
write_csv_atomic(de130, p("outputs", "de", "GSE130970_DESeq2_replication.csv"))
write_csv_atomic(de130[de130$significant, ], p("outputs", "de", "GSE130970_DESeq2_significant.csv"))
write_csv_atomic(de130_unadj, p("outputs", "de", "GSE130970_DESeq2_unadjusted_sensitivity.csv"))

saveRDS(list(counts = counts126, coldata = col126, vst = fit126$vst,
             de = de126, design = "~ group", population = "male-only",
             full_sensitivity = list(counts = counts126_all, coldata = col126_all,
                                     vst = fit126_full$vst, de = de126_full,
                                     design = "~ sex + group")),
        p("data", "derived", "GSE126848_rebuilt.rds"), compress = "xz")
saveRDS(list(counts = counts130, coldata = col130, vst = fit130$vst,
             de = de130, de_unadjusted = de130_unadj,
             design = "~ sex + age_z + group"),
        p("data", "derived", "GSE130970_rebuilt.rds"), compress = "xz")

common <- intersect(de126$Gene, de130$Gene)
a <- de126[match(common, de126$Gene), ]
b <- de130[match(common, de130$Gene), ]
valid <- is.finite(a$log2FoldChange) & is.finite(a$lfcSE) & a$lfcSE > 0 &
  is.finite(b$log2FoldChange) & is.finite(b$lfcSE) & b$lfcSE > 0
a <- a[valid, ]; b <- b[valid, ]; common <- common[valid]
w1 <- 1 / a$lfcSE^2
w2 <- 1 / b$lfcSE^2
meta_beta <- (w1 * a$log2FoldChange + w2 * b$log2FoldChange) / (w1 + w2)
meta_se <- sqrt(1 / (w1 + w2))
meta_z <- meta_beta / meta_se
meta_p <- 2 * pnorm(-abs(meta_z))
meta <- data.frame(
  Gene = common,
  log2FC_GSE126848 = a$log2FoldChange,
  SE_GSE126848 = a$lfcSE,
  log2FC_GSE130970 = b$log2FoldChange,
  SE_GSE130970 = b$lfcSE,
  meta_log2FC = meta_beta,
  meta_SE = meta_se,
  meta_z = meta_z,
  meta_p = meta_p,
  meta_fdr = p.adjust(meta_p, "BH"),
  same_direction = sign(a$log2FoldChange) == sign(b$log2FoldChange),
  stringsAsFactors = FALSE
)
meta <- meta[order(meta$meta_fdr, -abs(meta$meta_log2FC)), ]
write_csv_atomic(meta, p("outputs", "de", "MASLD_fixed_effect_meta_secondary.csv"))

benchmark_summary <- data.frame(
  metric = c(
    "n_genes_primary", "n_genes_replication", "n_common_finite",
    "n_significant_primary", "n_significant_replication",
    "all_gene_pearson_r", "all_gene_spearman_rho", "all_gene_sign_agreement",
    "primary_sig_replication_sign_agreement"
  ),
  value = c(
    nrow(de126), nrow(de130), length(common),
    sum(de126$significant), sum(de130$significant),
    cor(a$log2FoldChange, b$log2FoldChange, method = "pearson"),
    cor(a$log2FoldChange, b$log2FoldChange, method = "spearman"),
    mean(sign(a$log2FoldChange) == sign(b$log2FoldChange)),
    mean(sign(a$log2FoldChange[a$significant]) == sign(b$log2FoldChange[a$significant]))
  ),
  stringsAsFactors = FALSE
)
write_csv_atomic(benchmark_summary, p("outputs", "de", "independent_replication_summary.csv"))

plot_pca <- function(vst_mat, coldata, label) {
  pc <- prcomp(t(vst_mat), scale. = FALSE)
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  d <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], group = coldata$group,
                  sample = rownames(coldata))
  ggplot(d, aes(PC1, PC2, colour = group)) +
    geom_point(size = 2.6, alpha = 0.85) +
    scale_colour_manual(values = c(Control = "#0072B2", MASLD = "#D55E00")) +
    labs(title = label, x = sprintf("PC1 (%.1f%%)", ve[1]),
         y = sprintf("PC2 (%.1f%%)", ve[2]), colour = NULL) +
    theme_classic(base_size = 11)
}

for (item in list(
  list(plot = plot_pca(fit126$vst, col126, "GSE126848 rebuilt VST (male overlap population)"), name = "GSE126848_PCA"),
  list(plot = plot_pca(fit130$vst, col130, "GSE130970 rebuilt VST"), name = "GSE130970_PCA")
)) {
  ggsave(p("outputs", "de", paste0(item$name, ".png")), item$plot,
         width = 6.5, height = 5.2, dpi = 300)
  ggsave(p("outputs", "de", paste0(item$name, ".pdf")), item$plot,
         width = 6.5, height = 5.2)
}

de_contract <- list(
  created = format(Sys.time(), tz = "UTC", usetz = TRUE),
  spec_version = spec$spec_version,
  GSE126848 = list(
    n_control = sum(col126$group == "Control"),
    n_masld = sum(col126$group == "MASLD"),
    n_genes_after_filter = nrow(counts126),
    n_significant = sum(de126$significant),
    population = "male-only overlap population",
    design = "~ group",
    design_rank = fit126$design_rank,
    design_columns = fit126$design_columns,
    full_sensitivity_n = ncol(counts126_all),
    full_sensitivity_n_significant = sum(de126_full$significant),
    full_sensitivity_design = "~ sex + group"
  ),
  GSE130970 = list(
    n_control = sum(col130$group == "Control"),
    n_masld = sum(col130$group == "MASLD"),
    n_genes_after_filter = nrow(counts130),
    n_significant = sum(de130$significant),
    design = "~ sex + age_z + group",
    design_rank = fit130$design_rank,
    design_columns = fit130$design_columns,
    count_ingestion = "official tximport lengthScaledTPM counts rounded at DESeq2 boundary"
  ),
  independent_replication = as.list(setNames(benchmark_summary$value, benchmark_summary$metric))
)
write_json_atomic(de_contract, p("outputs", "de", "de_contract.json"))
save_session_info("01_masld_de")
cat("MASLD DE rebuild completed.\n")

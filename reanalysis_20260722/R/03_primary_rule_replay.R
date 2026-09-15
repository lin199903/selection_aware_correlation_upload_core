script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE)), "_common.R"))

required_packages(c("jsonlite", "digest"))

spec <- read_spec()
set.seed(spec$seed_primary)
ensure_dirs(p("outputs", "primary"), p("logs"))

fdr_cut <- as.numeric(spec$de_thresholds$fdr)
lfc_cut <- as.numeric(spec$de_thresholds$abs_log2fc)
n_perm <- as.integer(spec$rule_replay$n_permutations)
min_selected <- as.integer(spec$rule_replay$min_selected_genes)

# This list is fixed before inspecting the rebuilt results. It matches the
# exclusion policy used by the previous revision, but is now applied before
# all directional selection and all permutations.
y_chromosome_genes <- c(
  "UTY", "USP9Y", "DDX3Y", "KDM5D", "EIF1AY", "ZFY", "SRY",
  "NLGN4Y", "RPS4Y1", "RPS4Y2", "TSPY1", "RBMY1A1", "DAZ1",
  "PRKY", "AMELY", "TBL1Y", "PCDH11Y", "TMSB4Y", "VCY",
  "CDY1", "CDY2A", "HSFY1", "TXLNGY", "BPY2", "PRY"
)

input_paths <- c(
  MASLD_DE_primary = p("outputs", "de", "GSE126848_DESeq2_primary.csv"),
  MASLD_DE_full_sensitivity = p("outputs", "de", "GSE126848_DESeq2_full_sex_adjusted_sensitivity.csv"),
  MASLD_hubs_topology = p("outputs", "network", "MASLD_topology_hubs_primary.csv"),
  MASLD_hubs_phenotype_sensitivity = p("outputs", "network", "MASLD_phenotype_module_hubs_sensitivity.csv"),
  OSA_DE_primary = p("outputs", "osa_de", "GSE135917_limma_primary.csv"),
  OSA_DE_unadjusted_sensitivity = p("outputs", "osa_de", "GSE135917_limma_unadjusted_sensitivity.csv"),
  OSA_hubs_topology = p("outputs", "osa_network", "OSA_topology_hubs_primary.csv"),
  OSA_hubs_phenotype_sensitivity = p("outputs", "osa_network", "OSA_phenotype_module_hubs_sensitivity.csv")
)
if (any(!file.exists(input_paths))) {
  stop("Missing required inputs: ", paste(names(input_paths)[!file.exists(input_paths)], collapse = ", "))
}

read_masld_de <- function(path) {
  d <- read.csv(path, check.names = FALSE)
  required <- c("Gene", "log2FoldChange", "padj")
  if (!all(required %in% names(d))) stop("Unexpected MASLD DE schema: ", path)
  out <- data.frame(
    Gene = trimws(as.character(d$Gene)),
    MASLD_log2FC = as.numeric(d$log2FoldChange),
    MASLD_FDR = as.numeric(d$padj),
    stringsAsFactors = FALSE
  )
  out <- out[nzchar(out$Gene) & !duplicated(out$Gene), , drop = FALSE]
  out$MASLD_significant <- is.finite(out$MASLD_FDR) & out$MASLD_FDR < fdr_cut &
    is.finite(out$MASLD_log2FC) & abs(out$MASLD_log2FC) > lfc_cut
  out
}

read_hubs <- function(path) {
  d <- read.csv(path, check.names = FALSE)
  if (!"Gene" %in% names(d)) stop("Unexpected hub schema: ", path)
  unique(trimws(as.character(d$Gene[nzchar(trimws(as.character(d$Gene)))])))
}

read_osa_de <- function(path) {
osa_raw <- read.csv(path, check.names = FALSE)
if (!all(c("Gene", "logFC", "adj.P.Val") %in% names(osa_raw))) stop("Unexpected OSA DE schema")
out <- data.frame(
  Gene = trimws(as.character(osa_raw$Gene)),
  OSA_log2FC = as.numeric(osa_raw$logFC),
  OSA_FDR = as.numeric(osa_raw$adj.P.Val),
  stringsAsFactors = FALSE
)
out <- out[nzchar(out$Gene) & !duplicated(out$Gene), , drop = FALSE]
out$OSA_significant <- is.finite(out$OSA_FDR) & out$OSA_FDR < fdr_cut &
  is.finite(out$OSA_log2FC) & abs(out$OSA_log2FC) > lfc_cut
out
}
osa_primary <- read_osa_de(input_paths[["OSA_DE_primary"]])
osa_unadjusted <- read_osa_de(input_paths[["OSA_DE_unadjusted_sensitivity"]])
osa_hubs_topology <- read_hubs(input_paths[["OSA_hubs_topology"]])
osa_hubs_phenotype <- read_hubs(input_paths[["OSA_hubs_phenotype_sensitivity"]])

masld_primary <- read_masld_de(input_paths[["MASLD_DE_primary"]])
masld_full <- read_masld_de(input_paths[["MASLD_DE_full_sensitivity"]])
masld_hubs_topology <- read_hubs(input_paths[["MASLD_hubs_topology"]])
masld_hubs_phenotype <- read_hubs(input_paths[["MASLD_hubs_phenotype_sensitivity"]])

add_one_p <- function(null, observed, upper = TRUE) {
  valid <- is.finite(null)
  if (!any(valid) || !is.finite(observed)) return(NA_real_)
  exceed <- if (upper) null[valid] >= observed else abs(null[valid]) >= abs(observed)
  (1 + sum(exceed)) / (1 + sum(valid))
}

safe_cor <- function(x, y, method = "pearson") {
  if (length(x) < min_selected || length(y) != length(x) ||
      stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  suppressWarnings(as.numeric(stats::cor(x, y, method = method)))
}

build_pool <- function(masld_de, masld_hubs, osa_de, osa_hubs, analysis_id) {
  common <- merge(masld_de, osa_de, by = "Gene", all = FALSE, sort = FALSE)
  common <- common[is.finite(common$MASLD_log2FC) & is.finite(common$OSA_log2FC), , drop = FALSE]
  common <- common[!common$Gene %in% y_chromosome_genes, , drop = FALSE]

  common$OSA_hub_and_MASLD_DE <- common$Gene %in% osa_hubs & common$MASLD_significant
  common$MASLD_hub_and_OSA_DE <- common$Gene %in% masld_hubs & common$OSA_significant
  pool <- common[common$OSA_hub_and_MASLD_DE | common$MASLD_hub_and_OSA_DE, , drop = FALSE]
  pool$AnalysisID <- analysis_id
  pool$SameSignObserved <- pool$MASLD_log2FC * pool$OSA_log2FC > 0
  pool$EligibilityEvidence <- ifelse(
    pool$OSA_hub_and_MASLD_DE & pool$MASLD_hub_and_OSA_DE,
    "both_cross_qualification_routes",
    ifelse(pool$OSA_hub_and_MASLD_DE, "OSA_hub_x_MASLD_DE", "MASLD_hub_x_OSA_DE")
  )
  list(common = common, pool = pool)
}

run_rule_replay <- function(pool, seed, n_perm) {
  x <- pool$MASLD_log2FC
  y <- pool$OSA_log2FC
  observed_keep <- x * y > 0
  observed_n <- sum(observed_keep)
  if (nrow(pool) < min_selected || observed_n < min_selected) {
    stop("Too few genes for rule replay: pool=", nrow(pool), ", selected=", observed_n)
  }
  observed_pearson <- safe_cor(x[observed_keep], y[observed_keep], "pearson")
  observed_spearman <- safe_cor(x[observed_keep], y[observed_keep], "spearman")

  set.seed(seed)
  null_pearson <- rep(NA_real_, n_perm)
  null_spearman <- rep(NA_real_, n_perm)
  null_n <- integer(n_perm)
  for (i in seq_len(n_perm)) {
    yp <- sample(y, length(y), replace = FALSE)
    keep <- x * yp > 0
    null_n[i] <- sum(keep)
    if (null_n[i] >= min_selected) {
      null_pearson[i] <- safe_cor(x[keep], yp[keep], "pearson")
      null_spearman[i] <- safe_cor(x[keep], yp[keep], "spearman")
    }
  }

  exact <- which(null_n == observed_n & is.finite(null_pearson))
  tolerance <- 0L
  conditional <- exact
  while (length(conditional) < 100L && tolerance < nrow(pool)) {
    tolerance <- tolerance + 1L
    conditional <- which(abs(null_n - observed_n) <= tolerance & is.finite(null_pearson))
  }

  list(
    observed_keep = observed_keep,
    observed_n = observed_n,
    observed_pearson = observed_pearson,
    observed_spearman = observed_spearman,
    null_pearson = null_pearson,
    null_spearman = null_spearman,
    null_n = null_n,
    p_pearson_upper = add_one_p(null_pearson, observed_pearson, TRUE),
    p_spearman_upper = add_one_p(null_spearman, observed_spearman, TRUE),
    conditional_tolerance = tolerance,
    conditional_n = length(conditional),
    conditional_p_pearson_upper = add_one_p(null_pearson[conditional], observed_pearson, TRUE)
  )
}

run_naive_random_sets <- function(x, y, selected_n, n_perm, seed) {
  if (selected_n > length(x)) stop("Selected size exceeds comparator universe")
  set.seed(seed)
  out <- rep(NA_real_, n_perm)
  for (i in seq_len(n_perm)) {
    idx <- sample.int(length(x), selected_n, replace = FALSE)
    out[i] <- safe_cor(x[idx], y[idx], "pearson")
  }
  out
}

variants <- list(
  topology_primary = list(
    de = masld_primary,
    hubs = masld_hubs_topology,
    osa_de = osa_primary,
    osa_hubs = osa_hubs_topology,
    label = "Rebuilt covariate-adjusted OSA DE and topology-only hubs (primary)"
  ),
  masld_phenotype_hub_sensitivity = list(
    de = masld_primary,
    hubs = masld_hubs_phenotype,
    osa_de = osa_primary,
    osa_hubs = osa_hubs_topology,
    label = "MASLD phenotype-module hubs (sensitivity)"
  ),
  masld_full_sex_adjusted_de_sensitivity = list(
    de = masld_full,
    hubs = masld_hubs_topology,
    osa_de = osa_primary,
    osa_hubs = osa_hubs_topology,
    label = "MASLD full-cohort sex-adjusted DE (sensitivity)"
  ),
  osa_unadjusted_de_sensitivity = list(
    de = masld_primary,
    hubs = masld_hubs_topology,
    osa_de = osa_unadjusted,
    osa_hubs = osa_hubs_topology,
    label = "OSA unadjusted DE (sensitivity)"
  ),
  osa_phenotype_hub_sensitivity = list(
    de = masld_primary,
    hubs = masld_hubs_topology,
    osa_de = osa_primary,
    osa_hubs = osa_hubs_phenotype,
    label = "OSA adjusted phenotype-module hubs (sensitivity)"
  )
)

results <- list()
summary_rows <- list()
for (i in seq_along(variants)) {
  id <- names(variants)[i]
  v <- variants[[i]]
  built <- build_pool(v$de, v$hubs, v$osa_de, v$osa_hubs, id)
  replay <- run_rule_replay(built$pool, spec$seed_primary + i * 1000L, n_perm)
  whole_r <- safe_cor(built$common$MASLD_log2FC, built$common$OSA_log2FC, "pearson")

  selected <- built$pool[replay$observed_keep, , drop = FALSE]
  write_csv_atomic(built$pool, p("outputs", "primary", paste0(id, "_eligible_pool.csv")))
  write_csv_atomic(selected, p("outputs", "primary", paste0(id, "_same_sign_selected.csv")))

  results[[id]] <- list(
    label = v$label,
    common = built$common,
    pool = built$pool,
    selected = selected,
    replay = replay,
    whole_transcriptome_pearson = whole_r
  )
  summary_rows[[id]] <- data.frame(
    AnalysisID = id,
    Label = v$label,
    CommonGenes = nrow(built$common),
    EligiblePoolGenes = nrow(built$pool),
    SameSignSelectedGenes = replay$observed_n,
    WholeTranscriptomePearson = whole_r,
    ObservedPearson = replay$observed_pearson,
    RuleReplayNullMeanPearson = mean(replay$null_pearson, na.rm = TRUE),
    RuleReplayNullSDPearson = sd(replay$null_pearson, na.rm = TRUE),
    RuleReplayPUpperPearson = replay$p_pearson_upper,
    ObservedSpearman = replay$observed_spearman,
    RuleReplayNullMeanSpearman = mean(replay$null_spearman, na.rm = TRUE),
    RuleReplayPUpperSpearman = replay$p_spearman_upper,
    ValidPermutations = sum(is.finite(replay$null_pearson)),
    ConditionalSizeTolerance = replay$conditional_tolerance,
    ConditionalPermutations = replay$conditional_n,
    ConditionalPUpperPearson = replay$conditional_p_pearson_upper,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
rownames(summary_df) <- NULL

# Deliberately naive comparators are retained only to quantify how much the
# conclusion changes when the data-dependent same-sign rule is not replayed.
primary <- results[["topology_primary"]]
naive_common <- run_naive_random_sets(
  primary$common$MASLD_log2FC, primary$common$OSA_log2FC,
  nrow(primary$selected), n_perm, spec$seed_primary + 9001L
)
naive_pool <- run_naive_random_sets(
  primary$pool$MASLD_log2FC, primary$pool$OSA_log2FC,
  nrow(primary$selected), n_perm, spec$seed_primary + 9002L
)
naive_summary <- data.frame(
  Comparator = c("naive_random_paired_genes_from_common_transcriptome",
                 "naive_random_paired_genes_from_eligible_pool"),
  UniverseGenes = c(nrow(primary$common), nrow(primary$pool)),
  SetSize = nrow(primary$selected),
  ObservedPearson = primary$replay$observed_pearson,
  NullMeanPearson = c(mean(naive_common, na.rm = TRUE), mean(naive_pool, na.rm = TRUE)),
  NullSDPearson = c(sd(naive_common, na.rm = TRUE), sd(naive_pool, na.rm = TRUE)),
  PUpperAddOne = c(add_one_p(naive_common, primary$replay$observed_pearson, TRUE),
                   add_one_p(naive_pool, primary$replay$observed_pearson, TRUE)),
  ValidPermutations = c(sum(is.finite(naive_common)), sum(is.finite(naive_pool))),
  IsValidPrimaryTest = FALSE,
  stringsAsFactors = FALSE
)

root_norm <- normalizePath(ROOT, winslash = "/")
relative_to_root <- function(x) {
  value <- normalizePath(x, winslash = "/")
  prefix <- paste0(root_norm, "/")
  if (startsWith(value, prefix)) substring(value, nchar(prefix) + 1L) else value
}
manifest <- data.frame(
  Input = names(input_paths),
  RelativePath = vapply(input_paths, relative_to_root, character(1)),
  Bytes = as.numeric(file.info(input_paths)$size),
  SHA256 = vapply(input_paths, sha256_file, character(1)),
  ProvenanceStatus = rep("rebuilt", length(input_paths)),
  stringsAsFactors = FALSE
)

diag_dir <- p("outputs", "diagnostic_gene_label_replay")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
write_csv_atomic(summary_df, file.path(diag_dir, "gene_label_replay_summary.csv"))
write_csv_atomic(naive_summary, p("outputs", "primary", "naive_comparator_summary.csv"))
write_csv_atomic(manifest, file.path(diag_dir, "gene_label_replay_input_manifest.csv"))
saveRDS(list(
  spec_version = spec$spec_version,
  thresholds = list(fdr = fdr_cut, abs_log2fc = lfc_cut),
  y_chromosome_genes = y_chromosome_genes,
  variants = results,
  naive_comparators = list(common = naive_common, eligible_pool = naive_pool),
  manifest = manifest
), file.path(diag_dir, "gene_label_replay_null_distributions.rds"), compress = "xz")

plot_primary <- function(device, file) {
  if (device == "pdf") pdf(file, width = 10, height = 4.5)
  if (device == "png") png(file, width = 3000, height = 1350, res = 300)
  old <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  on.exit({par(old); dev.off()}, add = TRUE)
  hist(primary$replay$null_pearson, breaks = 60, col = "#D9EAF4", border = "white",
       xlab = "Pearson correlation", main = "Selection-aware rule replay")
  abline(v = primary$replay$observed_pearson, col = "#B2182B", lwd = 2, lty = 2)
  legend("topright", legend = sprintf("Observed r = %.3f\nP = %.4g",
                                      primary$replay$observed_pearson,
                                      primary$replay$p_pearson_upper),
         bty = "n", text.col = "#B2182B")
  hist(primary$replay$null_n, breaks = seq(min(primary$replay$null_n) - 0.5,
                                          max(primary$replay$null_n) + 0.5, by = 1),
       col = "#E8D9C5", border = "white", xlab = "Genes retained after same-sign rule",
       main = "Selected-set size under null")
  abline(v = primary$replay$observed_n, col = "#2166AC", lwd = 2, lty = 2)
}
plot_primary("pdf", file.path(diag_dir, "gene_label_replay.pdf"))
plot_primary("png", file.path(diag_dir, "gene_label_replay.png"))

contract <- list(
  created = format(Sys.time(), tz = "UTC", usetz = TRUE),
  spec_version = spec$spec_version,
  primary_analysis = "topology_primary",
  osa_inputs_rebuilt = TRUE,
  osa_input_status = "GSE135917 Study Group 1 DE and topology-only WGCNA rebuilt from the GPL6244 series matrix; Study Group 2 CPAP samples excluded from discovery.",
  y_exclusion_applied_before_selection_and_permutation = TRUE,
  rule_replayed_inside_every_permutation = TRUE,
  permutation_tail = spec$rule_replay$tail,
  add_one_correction = TRUE,
  n_permutations = n_perm,
  primary = as.list(summary_df[summary_df$AnalysisID == "topology_primary", , drop = FALSE]),
  naive_comparators_are_primary_tests = FALSE
)

stopifnot(
  !any(primary$pool$Gene %in% y_chromosome_genes),
  nrow(primary$selected) >= min_selected,
  sum(is.finite(primary$replay$null_pearson)) > 0.99 * n_perm,
  identical(spec$rule_replay$p_value_correction, "add_one"),
  identical(spec$rule_replay$selection_rule, "effect_a * effect_b > 0")
)
write_json_atomic(contract, file.path(diag_dir, "gene_label_replay_contract.json"))
save_session_info("03_primary_rule_replay")

print(summary_df)
print(naive_summary)
cat("Gene-label replay diagnostic and fixed-set comparator completed.\n")

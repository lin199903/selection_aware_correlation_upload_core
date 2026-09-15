script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE)), "_common.R"))

required_packages(c("WGCNA", "matrixStats", "jsonlite"))
suppressPackageStartupMessages({
  library(WGCNA)
  library(matrixStats)
})

spec <- read_spec()
set.seed(spec$seed_primary + 8L)
ensure_dirs(p("outputs", "osa_network"), p("logs"))

input_path <- p("data", "derived", "GSE135917_rebuilt.rds")
if (!file.exists(input_path)) stop("Run 08_osa_de.R first")
x <- readRDS(input_path)
expr <- x$discovery_expression
coldata <- x$discovery_coldata
stopifnot(identical(colnames(expr), rownames(coldata)), ncol(expr) == 18L)

mad_values <- rowMads(expr)
names(mad_values) <- rownames(expr)
mad_values <- mad_values[is.finite(mad_values) & mad_values > 0]
n_keep <- min(as.integer(spec$network$top_mad_genes), length(mad_values))
network_genes <- names(sort(mad_values, decreasing = TRUE))[seq_len(n_keep)]
datExpr <- t(expr[network_genes, , drop = FALSE])
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
stopifnot(nrow(datExpr) == 18L, ncol(datExpr) >= 1000L, !anyDuplicated(colnames(datExpr)))

allowWGCNAThreads(nThreads = max(1L, min(8L, parallel::detectCores() - 1L)))
powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(
  datExpr, powerVector = powers, networkType = spec$network$network_type,
  corFnc = "bicor", verbose = 0
)
fit <- sft$fitIndices
signed_r2 <- -sign(fit[, 3]) * fit[, 2]
eligible_power <- fit[, 1][is.finite(signed_r2) & signed_r2 >= 0.80]
power <- if (length(eligible_power)) min(eligible_power) else {
  candidate <- fit[, 1][which.max(ifelse(is.finite(signed_r2), signed_r2, -Inf))]
  if (!is.finite(candidate)) 6 else candidate
}

plot_sft <- function(device, file) {
  if (device == "pdf") pdf(file, width = 10, height = 4.5)
  if (device == "png") png(file, width = 3000, height = 1350, res = 300)
  old <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  on.exit({par(old); dev.off()}, add = TRUE)
  plot(fit[, 1], signed_r2, type = "n", xlab = "Soft-threshold power",
       ylab = "Signed scale-free fit (R²)", main = "Scale independence")
  text(fit[, 1], signed_r2, labels = fit[, 1], col = "#B2182B")
  abline(h = 0.80, lty = 2, col = "grey40")
  abline(v = power, lty = 3, col = "#2166AC")
  plot(fit[, 1], fit[, 5], type = "n", xlab = "Soft-threshold power",
       ylab = "Mean connectivity", main = "Mean connectivity")
  text(fit[, 1], fit[, 5], labels = fit[, 1], col = "#2166AC")
  abline(v = power, lty = 3, col = "#2166AC")
}
plot_sft("pdf", p("outputs", "osa_network", "soft_threshold.pdf"))
plot_sft("png", p("outputs", "osa_network", "soft_threshold.png"))

net <- blockwiseModules(
  datExpr,
  power = power,
  networkType = spec$network$network_type,
  TOMType = spec$network$tom_type,
  corType = "bicor",
  maxBlockSize = ncol(datExpr) + 1L,
  minModuleSize = as.integer(spec$network$min_module_size),
  mergeCutHeight = spec$network$merge_cut_height,
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  reassignThreshold = 0,
  saveTOMs = FALSE,
  verbose = 1
)
module_colors <- net$colors
names(module_colors) <- colnames(datExpr)
MEs <- orderMEs(moduleEigengenes(datExpr, module_colors)$eigengenes)
kme <- signedKME(datExpr, MEs, outputColumnName = "kME", corFnc = "bicor")

module_table <- data.frame(
  Gene = colnames(datExpr),
  Module = unname(module_colors[colnames(datExpr)]),
  MAD = mad_values[colnames(datExpr)],
  stringsAsFactors = FALSE
)
module_table$kME <- vapply(seq_len(nrow(module_table)), function(i) {
  col <- paste0("kME", module_table$Module[i])
  if (col %in% colnames(kme)) kme[i, col] else NA_real_
}, numeric(1))

hub_rows <- lapply(setdiff(unique(module_table$Module), "grey"), function(mod) {
  d <- module_table[module_table$Module == mod &
                      is.finite(module_table$kME) &
                      abs(module_table$kME) >= spec$network$hub_abs_kme, , drop = FALSE]
  d <- d[order(-abs(d$kME), -d$MAD), , drop = FALSE]
  head(d, as.integer(spec$network$hub_max_per_module))
})
topology_hubs <- do.call(rbind, hub_rows)
if (is.null(topology_hubs)) topology_hubs <- module_table[0, ]
topology_hubs$HubRule <- "topology_only_abs_kME"

trait <- as.numeric(coldata[rownames(datExpr), "group"] == "OSA")
module_cor <- bicor(MEs, trait, use = "pairwise.complete.obs")
module_p <- corPvalueStudent(module_cor, nrow(datExpr))
module_adjusted <- do.call(rbind, lapply(seq_len(ncol(MEs)), function(j) {
  d <- data.frame(
    ME = MEs[, j],
    age_z = coldata[rownames(MEs), "age_z"],
    bmi_z = coldata[rownames(MEs), "bmi_z"],
    sex = coldata[rownames(MEs), "sex"],
    group = coldata[rownames(MEs), "group"]
  )
  fit_lm <- lm(ME ~ age_z + bmi_z + sex + group, data = d)
  co <- summary(fit_lm)$coefficients["groupOSA", ]
  data.frame(
    ModuleEigengene = colnames(MEs)[j], Module = sub("^ME", "", colnames(MEs)[j]),
    AdjustedGroupBeta = unname(co["Estimate"]),
    AdjustedGroupSE = unname(co["Std. Error"]),
    AdjustedGroupP = unname(co["Pr(>|t|)"]), stringsAsFactors = FALSE
  )
}))
module_trait <- merge(
  data.frame(
    ModuleEigengene = rownames(module_cor),
    Module = sub("^ME", "", rownames(module_cor)),
    RawCorrelation = as.numeric(module_cor[, 1]),
    RawPvalue = as.numeric(module_p[, 1]), stringsAsFactors = FALSE
  ),
  module_adjusted, by = c("ModuleEigengene", "Module"), sort = FALSE
)
module_trait$AdjustedFDR <- p.adjust(module_trait$AdjustedGroupP, "BH")
sig_modules <- module_trait$Module[module_trait$AdjustedGroupP < 0.05 & module_trait$Module != "grey"]
phenotype_hubs <- topology_hubs[topology_hubs$Module %in% sig_modules, , drop = FALSE]
phenotype_hubs$HubRule <- "adjusted_phenotype_module_nominal_p_lt_0.05_sensitivity"

write_csv_atomic(module_table, p("outputs", "osa_network", "OSA_module_assignments.csv"))
write_csv_atomic(topology_hubs, p("outputs", "osa_network", "OSA_topology_hubs_primary.csv"))
write_csv_atomic(phenotype_hubs, p("outputs", "osa_network", "OSA_phenotype_module_hubs_sensitivity.csv"))
write_csv_atomic(module_trait, p("outputs", "osa_network", "OSA_module_trait_sensitivity.csv"))

plot_dendro <- function(device, file) {
  if (device == "pdf") pdf(file, width = 12, height = 6)
  if (device == "png") png(file, width = 3600, height = 1800, res = 300)
  on.exit(dev.off(), add = TRUE)
  plotDendroAndColors(
    net$dendrograms[[1]], module_colors[net$blockGenes[[1]]], "Modules",
    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
    main = "GSE135917 Study Group 1 topology-only network"
  )
}
plot_dendro("pdf", p("outputs", "osa_network", "module_dendrogram.pdf"))
plot_dendro("png", p("outputs", "osa_network", "module_dendrogram.png"))

saveRDS(list(
  spec_version = spec$spec_version,
  population = "GSE135917 Study Group 1: 8 controls and 10 OSA cases",
  power = power,
  fit_indices = fit,
  module_colors = module_colors,
  module_eigengenes = MEs,
  module_table = module_table,
  topology_hubs = topology_hubs,
  phenotype_hubs = phenotype_hubs,
  module_trait = module_trait,
  network_parameters = spec$network
), p("outputs", "osa_network", "OSA_network_rebuilt.rds"), compress = "xz")

contract <- list(
  created = format(Sys.time(), tz = "UTC", usetz = TRUE),
  spec_version = spec$spec_version,
  population = "GSE135917 Study Group 1",
  n_samples = nrow(datExpr),
  n_input_genes = ncol(datExpr),
  soft_power = power,
  n_modules_including_grey = length(unique(module_colors)),
  n_topology_hubs = nrow(topology_hubs),
  n_adjusted_phenotype_module_hubs_sensitivity = nrow(phenotype_hubs),
  duplicated_gene_ids = anyDuplicated(module_table$Gene),
  hub_definition = "all non-grey modules; abs(kME)>=0.80; top 20 per module",
  phenotype_module_filter_used_in_primary = FALSE
)
stopifnot(contract$n_samples == 18L,
          contract$n_input_genes <= spec$network$top_mad_genes,
          contract$duplicated_gene_ids == 0L,
          nrow(topology_hubs) > 0L)
write_json_atomic(contract, p("outputs", "osa_network", "osa_network_contract.json"))
save_session_info("09_osa_network")
cat("OSA topology-only WGCNA rebuild completed.\n")

# 12_osa_recovery_benchmark.R
# Null-calibrated calibration-recovery benchmark for the disease-label replay.
# The z-score machinery (correlated_z, signal_index, selection_mask,
# selected_cor, block and gene-label permutations, add-one Monte Carlo P)
# follows the v107 factorial screen (04_systematic_simulation.R) so that the
# calibration cells are directly comparable. New shared-signal generators are
# added: a single-arm null (MASLD-like signal, no OSA-side signal), sparse
# shared spike-ins at stronger effect sizes, dense proportional shared effects
# (beta_B = c * theta + noise), and directionally heterogeneous shared effects.

script_arg <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE)), "_common.R"))

required_packages(c("jsonlite", "matrixStats"))

spec <- read_spec()
rb <- spec$benchmarks$recovery_benchmark
alpha <- as.numeric(spec$benchmarks$alpha)
seed_base <- as.integer(spec$seed_benchmark) + as.integer(rb$seed_offset)
outer_R <- as.integer(rb$replicates)
inner_B <- as.integer(rb$reference_replicates)
n_group_anchor <- as.integer(rb$group_n_anchor)
rho_anchor <- as.numeric(rb$block_rho_anchor)
confounding_levels <- as.character(rb$confounding_levels)

smoke <- "--smoke" %in% commandArgs(trailingOnly = TRUE)
if (smoke) {
  outer_R <- 20L
  inner_B <- 99L
}

G_screen <- 120L
block_size_screen <- 10L
selection_fraction <- 0.50

out_dir <- p("outputs", "recovery_benchmark")
ensure_dirs(out_dir, p("logs"))

block_ids <- function(G, block_size) {
  rep(seq_len(ceiling(G / block_size)), each = block_size)[seq_len(G)]
}

correlated_z <- function(G, R, rho, block_size) {
  bid <- block_ids(G, block_size)
  B <- max(bid)
  block_component <- matrix(stats::rnorm(B * R), nrow = B, ncol = R)
  gene_component <- matrix(stats::rnorm(G * R), nrow = G, ncol = R)
  sqrt(rho) * block_component[bid, , drop = FALSE] + sqrt(1 - rho) * gene_component
}

signal_index <- function(G, sparsity, block_size) {
  k <- as.integer(round(G * sparsity))
  if (k <= 0L) return(integer())
  bid <- block_ids(G, block_size)
  ord <- order(ave(seq_len(G), bid, FUN = seq_along), bid)
  sort(ord[seq_len(k)])
}

confounder_loading_index <- function(G, block_size) {
  bid <- block_ids(G, block_size)
  within <- ave(seq_len(G), bid, FUN = seq_along)
  which(within <= max(1L, floor(block_size / 4L)))
}

top_fraction <- function(tmat, fraction = selection_fraction) {
  ranks <- matrixStats::colRanks(abs(tmat), ties.method = "first", preserveShape = TRUE)
  ranks > nrow(tmat) * (1 - fraction)
}

selection_mask <- function(tA, tB) {
  top_fraction(tA) & top_fraction(tB) & (sign(tA) == sign(tB))
}

selected_cor <- function(betaA, betaB, selected) {
  n <- colSums(selected)
  sA <- colSums(betaA * selected)
  sB <- colSums(betaB * selected)
  sAA <- colSums(betaA * betaA * selected)
  sBB <- colSums(betaB * betaB * selected)
  sAB <- colSums(betaA * betaB * selected)
  covAB <- sAB - sA * sB / pmax(n, 1)
  varA <- sAA - sA * sA / pmax(n, 1)
  varB <- sBB - sB * sB / pmax(n, 1)
  ans <- covAB / sqrt(pmax(varA * varB, 0))
  ans[n < 4L | !is.finite(ans)] <- 0
  pmax(-1, pmin(1, ans))
}

permute_genes_by_column <- function(m) {
  out <- m
  for (j in seq_len(ncol(m))) out[, j] <- m[sample.int(nrow(m)), j]
  out
}

permute_blocks_by_column <- function(m, block_size) {
  G <- nrow(m)
  bid <- block_ids(G, block_size)
  B <- max(bid)
  out <- m
  for (j in seq_len(ncol(m))) {
    new_order <- sample.int(B)
    idx <- unlist(lapply(new_order, function(b) which(bid == b)), use.names = FALSE)
    out[, j] <- m[idx, j]
  }
  out
}

mc_p <- function(observed, reference) {
  (1 + vapply(observed, function(x) sum(reference >= x), integer(1))) / (length(reference) + 1)
}

wilson_ci <- function(x, n, level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  phat <- x / n
  den <- 1 + z^2 / n
  centre <- (phat + z^2 / (2 * n)) / den
  half <- z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2)) / den
  c(lower = max(0, centre - half), upper = min(1, centre + half))
}

simulate_effects <- function(R, n_group, rho, sparsity, effect_sd, confounding,
                             G = G_screen, block_size = block_size_screen) {
  se <- sqrt(2 / n_group)
  z <- correlated_z(G, R, rho, block_size)
  beta <- se * z
  idx <- signal_index(G, sparsity, block_size)
  if (length(idx) && effect_sd > 0) beta[idx, ] <- beta[idx, , drop = FALSE] + effect_sd
  if (identical(confounding, "misspecified")) {
    cidx <- confounder_loading_index(G, block_size)
    loadings <- seq(-0.60, 0.60, length.out = length(cidx))
    beta[cidx, ] <- beta[cidx, , drop = FALSE] + loadings
  }
  list(beta = beta, t = beta / se)
}

generate_scenario <- function(type, R, n_group, rho, confounding, seed,
                              sparsity = NULL, effect_sd = NULL,
                              c = NULL, h = NULL, flip_fraction = NULL,
                              single_arm = FALSE) {
  se <- sqrt(2 / n_group)
  set.seed(seed)
  theta <- if (type == "shared") correlated_z(G_screen, R, rho, block_size_screen) else NULL
  eta <- NULL
  if (type == "shared" && !is.null(flip_fraction) && flip_fraction > 0) {
    set.seed(seed + 1L)
    eta <- (stats::runif(G_screen) > flip_fraction) * 2 - 1
  }
  make_arm <- function(which) {
    z <- correlated_z(G_screen, R, rho, block_size_screen)
    beta <- se * z
    if (!is.null(theta)) {
      sign_mult <- if (which == "B") (if (is.null(eta)) c else c * eta) else 1
      beta <- beta + h * se * theta * sign_mult
    }
    idx <- if (type == "spike") signal_index(G_screen, sparsity, block_size_screen) else integer()
    apply_signal <- !single_arm || which == "A"
    if (length(idx) && !is.null(effect_sd) && effect_sd > 0 && apply_signal) {
      beta[idx, ] <- beta[idx, , drop = FALSE] + effect_sd
    }
    if (identical(confounding, "misspecified") && apply_signal) {
      cidx <- confounder_loading_index(G_screen, block_size_screen)
      loadings <- seq(-0.60, 0.60, length.out = length(cidx))
      beta[cidx, ] <- beta[cidx, , drop = FALSE] + loadings
    }
    list(beta = beta, t = beta / se)
  }
  list(A1 = make_arm("A"), B1 = make_arm("B"), A2 = make_arm("A"), B2 = make_arm("B"))
}

make_references <- function(A, B, n_group, rho, confounding) {
  selected_observed <- selection_mask(A$t, B$t)
  B_gene_beta <- permute_genes_by_column(B$beta)
  B_gene_t <- B_gene_beta / sqrt(2 / n_group)
  naive <- selected_cor(A$beta, B_gene_beta, selected_observed)
  replay <- selected_cor(A$beta, B_gene_beta, selection_mask(A$t, B_gene_t))
  B_block_beta <- permute_blocks_by_column(B$beta, block_size_screen)
  B_block_t <- B_block_beta / sqrt(2 / n_group)
  block <- selected_cor(A$beta, B_block_beta, selection_mask(A$t, B_block_t))
  A_null <- simulate_effects(inner_B, n_group, rho, 0, 0, confounding)
  B_null <- simulate_effects(inner_B, n_group, rho, 0, 0, confounding)
  full_two_arm <- selected_cor(A_null$beta, B_null$beta, selection_mask(A_null$t, B_null$t))
  list(
    naive_fixed_set_gene_label = naive,
    gene_label_complete_rule_replay = replay,
    block_preserving_rule_replay = block,
    full_two_arm_synthetic_replay = full_two_arm
  )
}

run_observed <- function(A, B, refs, R) {
  selected <- selection_mask(A$t, B$t)
  observed <- selected_cor(A$beta, B$beta, selected)
  pvals <- lapply(refs, function(x) mc_p(observed, x))
  list(
    observed = observed,
    selected_n = colSums(selected),
    p = pvals,
    reject = lapply(pvals, function(x) x <= alpha)
  )
}

run_sample_splitting <- function(A1, B1, A2, B2, R) {
  selected1 <- selection_mask(A1$t, B1$t)
  n_sel <- colSums(selected1)
  r2 <- selected_cor(A2$beta, B2$beta, selected1)
  zf <- atanh(pmax(-0.9999, pmin(0.9999, r2))) * sqrt(pmax(n_sel - 3L, 0))
  p <- stats::pnorm(zf, lower.tail = FALSE)
  list(observed = r2, selected_n = n_sel, p = p, reject = p <= alpha)
}

scenarios <- list(
  list(mechanism = "global_null", type = "spike", sparsity = 0, effect_sd = 0, single_arm = FALSE),
  list(mechanism = "single_arm_null", type = "spike", sparsity = 0.1, effect_sd = 0.8, single_arm = TRUE),
  list(mechanism = "sparse_shared", type = "spike", sparsity = 0.1, effect_sd = 0.8, single_arm = FALSE),
  list(mechanism = "sparse_shared", type = "spike", sparsity = 0.1, effect_sd = 1.2, single_arm = FALSE),
  list(mechanism = "proportional_shared", type = "shared", c = 0.5, h = 0.75, single_arm = FALSE),
  list(mechanism = "proportional_shared", type = "shared", c = 0.5, h = 1.5, single_arm = FALSE),
  list(mechanism = "proportional_shared", type = "shared", c = 0.5, h = 2.5, single_arm = FALSE),
  list(mechanism = "proportional_shared", type = "shared", c = 1.0, h = 0.75, single_arm = FALSE),
  list(mechanism = "proportional_shared", type = "shared", c = 1.0, h = 1.5, single_arm = FALSE),
  list(mechanism = "proportional_shared", type = "shared", c = 1.0, h = 2.5, single_arm = FALSE),
  list(mechanism = "heterogeneous", type = "shared", c = 1.0, h = 1.5, flip_fraction = 0.25, single_arm = FALSE),
  list(mechanism = "heterogeneous", type = "shared", c = 1.0, h = 1.5, flip_fraction = 0.5, single_arm = FALSE)
)

method_names <- c(
  "naive_fixed_set_gene_label",
  "gene_label_complete_rule_replay",
  "block_preserving_rule_replay",
  "full_two_arm_synthetic_replay",
  "sample_splitting_fixed_set"
)

scenario_id <- function(s, ci) {
  base <- s$mechanism
  if (base == "proportional_shared") {
    paste0("prop_shared_c", s$c, "_h", s$h)
  } else if (base == "heterogeneous") {
    paste0("hetero_flip", s$flip_fraction)
  } else if (base == "sparse_shared") {
    paste0("sparse_shared_", s$effect_sd)
  } else {
    base
  }
}

start_time <- Sys.time()
summary_rows <- list()
ref_diag_rows <- list()
raw_rows <- list()

for (ci in seq_along(confounding_levels)) {
  conf <- confounding_levels[ci]
  for (si in seq_along(scenarios)) {
    s <- scenarios[[si]]
    id <- scenario_id(s, ci)
    seed_s <- seed_base + (si - 1L) * 100L + (ci - 1L) * 10L + 1L
    message(sprintf("[%s] n=%d rho=%.1f conf=%s rep=%d ref=%d",
                    id, n_group_anchor, rho_anchor, conf, outer_R, inner_B))
    dat <- generate_scenario(
      type = s$type, R = outer_R, n_group = n_group_anchor, rho = rho_anchor,
      confounding = conf, seed = seed_s,
      sparsity = s$sparsity, effect_sd = s$effect_sd,
      c = s$c, h = s$h, flip_fraction = s$flip_fraction,
      single_arm = s$single_arm
    )
    refs <- make_references(dat$A1, dat$B1, n_group_anchor, rho_anchor, confounding = conf)
    obs <- run_observed(dat$A1, dat$B1, refs, outer_R)
    split <- run_sample_splitting(dat$A1, dat$B1, dat$A2, dat$B2, outer_R)

    ref_q95 <- vapply(refs, function(x) unname(stats::quantile(x, 0.95)), numeric(1))
    for (m in names(refs)) {
      x <- sum(obs$reject[[m]])
      ci_w <- wilson_ci(x, outer_R)
      ref_n <- if (m == "full_two_arm_synthetic_replay") inner_B else outer_R
      ref_construction <- if (m == "full_two_arm_synthetic_replay") {
        "independent_synthetic_null_draws"
      } else {
        "pooled_cross_replicate_permutation"
      }
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        mechanism = s$mechanism, scenario = id, confounding = conf,
        c = ifelse(is.null(s$c), NA_real_, s$c),
        h = ifelse(is.null(s$h), NA_real_, s$h),
        sparsity = ifelse(is.null(s$sparsity), NA_real_, s$sparsity),
        effect_sd = ifelse(is.null(s$effect_sd), NA_real_, s$effect_sd),
        flip_fraction = ifelse(is.null(s$flip_fraction), NA_real_, s$flip_fraction),
        replicates = outer_R, reference_replicates = ref_n,
        reference_construction = ref_construction,
        method = m, rejection_rate = x / outer_R,
        wilson_lower = ci_w[["lower"]], wilson_upper = ci_w[["upper"]],
        mcse = sqrt((x / outer_R) * (1 - x / outer_R) / outer_R),
        covers_nominal_0.05 = ci_w[["lower"]] <= alpha && ci_w[["upper"]] >= alpha,
        mean_selected_n = mean(obs$selected_n),
        median_statistic = stats::median(obs$observed),
        reference_q95 = ref_q95[[m]],
        separation = stats::median(obs$observed) - ref_q95[[m]],
        stringsAsFactors = FALSE
      )
      ref_diag_rows[[length(ref_diag_rows) + 1L]] <- data.frame(
        mechanism = s$mechanism, scenario = id, confounding = conf,
        method = m, reference_median = stats::median(refs[[m]]),
        reference_q95 = ref_q95[[m]],
        stringsAsFactors = FALSE
      )
      raw_rows[[length(raw_rows) + 1L]] <- data.frame(
        mechanism = s$mechanism, scenario = id, confounding = conf,
        replicate = seq_len(outer_R), method = m,
        selected_n = obs$selected_n, statistic = obs$observed,
        p_value = obs$p[[m]], reject = obs$reject[[m]],
        stringsAsFactors = FALSE
      )
    }
    x <- sum(split$reject)
    ci_w <- wilson_ci(x, outer_R)
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      mechanism = s$mechanism, scenario = id, confounding = conf,
      c = ifelse(is.null(s$c), NA_real_, s$c),
      h = ifelse(is.null(s$h), NA_real_, s$h),
      sparsity = ifelse(is.null(s$sparsity), NA_real_, s$sparsity),
      effect_sd = ifelse(is.null(s$effect_sd), NA_real_, s$effect_sd),
      flip_fraction = ifelse(is.null(s$flip_fraction), NA_real_, s$flip_fraction),
      replicates = outer_R, reference_replicates = NA_integer_,
      reference_construction = "no_reference_fisher_z",
      method = "sample_splitting_fixed_set", rejection_rate = x / outer_R,
      wilson_lower = ci_w[["lower"]], wilson_upper = ci_w[["upper"]],
      mcse = sqrt((x / outer_R) * (1 - x / outer_R) / outer_R),
      covers_nominal_0.05 = ci_w[["lower"]] <= alpha && ci_w[["upper"]] >= alpha,
      mean_selected_n = mean(split$selected_n),
      median_statistic = stats::median(split$observed),
      reference_q95 = NA_real_,
      separation = NA_real_,
      stringsAsFactors = FALSE
    )
    raw_rows[[length(raw_rows) + 1L]] <- data.frame(
      mechanism = s$mechanism, scenario = id, confounding = conf,
      replicate = seq_len(outer_R), method = "sample_splitting_fixed_set",
      selected_n = split$selected_n, statistic = split$observed,
      p_value = split$p, reject = split$reject,
      stringsAsFactors = FALSE
    )
  }
}

summary_df <- do.call(rbind, summary_rows)
ref_diag_df <- do.call(rbind, ref_diag_rows)
raw_df <- do.call(rbind, raw_rows)
row.names(summary_df) <- NULL
row.names(ref_diag_df) <- NULL
row.names(raw_df) <- NULL

write_csv_atomic(summary_df, file.path(out_dir, "recovery_benchmark_summary.csv"))
write_csv_atomic(ref_diag_df, file.path(out_dir, "recovery_benchmark_reference_diagnostics.csv"))
saveRDS(raw_df, file.path(out_dir, "recovery_benchmark_raw_replicates.rds"), compress = "xz")

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
contract <- list(
  spec_version = spec$spec_version,
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  seed_base = seed_base,
  alpha = alpha,
  G_screen = G_screen,
  block_size = block_size_screen,
  selection_fraction = selection_fraction,
  group_n_anchor = n_group_anchor,
  block_rho_anchor = rho_anchor,
  confounding_levels = confounding_levels,
  replicates = outer_R,
  reference_replicates = inner_B,
  smoke_test = smoke,
  scenario_count = length(scenarios) * length(confounding_levels),
  elapsed_minutes = round(elapsed, 2),
  outputs = list(
    summary = sha256_file(file.path(out_dir, "recovery_benchmark_summary.csv")),
    reference_diagnostics = sha256_file(file.path(out_dir, "recovery_benchmark_reference_diagnostics.csv")),
    raw_replicates = sha256_file(file.path(out_dir, "recovery_benchmark_raw_replicates.rds"))
  )
)
write_json_atomic(contract, file.path(out_dir, "recovery_benchmark_contract.json"))
save_session_info("12_osa_recovery_benchmark")

message(sprintf("Recovery benchmark complete: %d scenarios, %.1f minutes elapsed.",
                length(scenarios) * length(confounding_levels), elapsed))

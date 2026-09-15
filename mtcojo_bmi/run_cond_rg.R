# Conditional genetic correlation after mtCOJO adjustment (BMI-conditioning)
# Reuses the exact main-text pipeline: pleioh2g::ldsc_rg, 200 jackknife blocks

suppressPackageStartupMessages({
  library(pleioh2g)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript run_cond_rg.R <workspace>\n  workspace: directory containing munged/ and ld_ref/ subdirectories")
}
workspace <- normalizePath(args[[1]])
sumstats_dir <- file.path(workspace, "munged")
ld_path  <- file.path(workspace, "ld_ref", "LDscore")
wld_path <- file.path(workspace, "ld_ref", "1000G_Phase3_weights_hm3_no_MHC")

run_pair <- function(run_tag, t1, t2, prev1, prev2) {
  cat("\n========================================\n")
  cat("Conditional rg:", run_tag, "\n")
  munged <- setNames(lapply(c(t1, t2), function(t) {
    fread(file.path(sumstats_dir, paste0(t, ".sumstats.gz")))
  }), paste0("GWAS_", c(t1, t2)))
  sample_prev <- c(prev1[1], prev2[1])
  population_prev <- c(prev1[2], prev2[2])
  res <- tryCatch({
    pleioh2g::ldsc_rg(
      munged_sumstats = munged,
      sample_prev = sample_prev,
      population_prev = population_prev,
      ld = ld_path,
      wld = wld_path,
      n_blocks = 200,
      chisq_max = NA,
      chr_filter = 1:22
    )
  }, error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) return(invisible(NULL))
  rg <- as.data.frame(res$rg)
  print(rg)
  out <- file.path(workspace, "mtcojo_bmi", paste0("cond_rg_", run_tag, ".csv"))
  write.csv(rg, out, row.names = FALSE)
  cat("Saved:", out, "\n")
  saveRDS(res, file.path(workspace, "mtcojo_bmi", paste0("cond_rg_", run_tag, ".rds")))
}

# Prevalence: NAFLD sample 4761/377762, pop 0.25; OSA sample 16761/202382, pop 0.10 (main-text contract)
NAFLD_P <- c(4761 / 377762, 0.25)
OSA_P   <- c(16761 / 202382, 0.10)
BMI_P   <- c(NA, NA)

# A. Main test: OSA conditional on BMI  vs  NAFLD
run_pair("OSA_cond_BMI_x_NAFLD", "OSA_cond_BMI", "NAFLD", BMI_P, NAFLD_P)
# B. Symmetric: NAFLD conditional on OSA+BMI  vs  OSA
run_pair("NAFLD_cond_OSA_BMI_x_OSA", "NAFLD_cond_OSA_BMI", "OSA", c(NAFLD_P[1], NAFLD_P[2]), OSA_P)
# C. Reproduction control: unadjusted OSA vs NAFLD (should reproduce r_g = 0.442)
run_pair("OSA_x_NAFLD_repro", "OSA", "NAFLD", OSA_P, NAFLD_P)

cat("\n=== All conditional rg analyses done ===\n")

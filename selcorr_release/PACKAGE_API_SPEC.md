# PACKAGE_API_SPEC — selcorr 0.1

## 1. Scope

`selcorr` implements one core object: a construction-aware reference for a
correlation computed after agreement-based feature selection. Two further
functions are subordinate: `diagnose()` (recovery boundary) and `validate()`
(portability in independent data). No composite score, no renamed statistic and
no additional framework are defined in 0.1.

## 2. Primary API

```r
fit <- selcorr(x, y, selection = selection_spec("de_anchor"), method = "replay", B = 2000)
summary(fit)
plot(fit)
```

### `selcorr(x, y, selection, anchor, method, B, seed, checkpoint, checkpoint_every, resume, parallel, keep_data, verbose)`

| Argument | Meaning |
|---|---|
| `x`, `y` | Arms. Effect-level: data frame with `gene`, `effect`, optional `fdr`, `hub`. Sample-level: list with `expr`, `group`, optional `covariates`, `design_full`, `design_nuisance`, `estimator` (`"limma"`), `fit_fun`. |
| `selection` | `selection_spec()` object. Required and recorded in `fit$spec`. |
| `anchor` | `"x"` (default) or `"y"` — which arm stays fixed. |
| `method` | `"replay"` (sample-level, primary), `"gene_label"` (effect-level), `"fixed_set"` (diagnostic comparator only). |
| `B` | Reference draws. `B = 0` returns the observed statistics without a reference. |
| `seed` | Integer. Draw sequences are generated before evaluation, so identical seeds reproduce identical draws in serial and parallel runs. |
| `checkpoint`, `checkpoint_every`, `resume` | Long-run persistence. Checkpoints store iteration, converged flag, statistic, eligible/selected sizes, common-gene count and error text. |
| `parallel` | Use `future.apply` with a multisession plan. |
| `keep_data` | Retain the inputs needed by `diagnose()` and `validate()`. |

Returned fields: `observed_r`, `n_common`, `n_eligible`, `n_selected`,
`selected_fraction`, `null_mean`, `null_sd`, `null_median`, `null_q025`,
`null_q95`, `critical_value`, `observed_percentile`, `z_position`, `p_value`,
`monte_carlo_se`, `n_reference_valid`, `n_reference_failed`, `selected_genes`,
`replay_statistics`, `spec`, `seed`, `B`, `method`, `diagnostic_only`.

### `selection_spec(eligibility, fdr_cut, abs_effect_cut, direction, min_selected, exclude)`

| Argument | Meaning |
|---|---|
| `eligibility` | `"de_anchor"` (pool = anchored arm's DE genes), `"hub_de"` (`(hub_replayed AND DE_anchored) OR (hub_anchored AND DE_replayed)`), or a function of the paired table. |
| `fdr_cut`, `abs_effect_cut` | Magnitude rule: `FDR < fdr_cut` and `|effect| > abs_effect_cut`. Defaults 0.05 and 0.3, matching the archived analyses. |
| `direction` | Only `"same_sign"` in 0.1. |
| `min_selected` | Minimum selected-set size for a defined correlation (default 4). |
| `exclude` | Genes removed by convention (Y chromosome by default). |

### `summary(fit)`

Prints the observed correlation, selected fraction, reference type, reference
mean/SD/95% critical value, observed percentile, upper-tail p value and z
position, an interpretation line restricted to "excess alignment relative to
this workflow-specific reference: detected / not detected", and a boundary
paragraph stating that the result concerns the selected-correlation endpoint
under the specified workflow.

### `plot(fit)`

One figure: the reference distribution with vertical lines for the observed
correlation, the reference centre and the 95% critical value. A diagnostic
fixed-set fit is labelled as such.

## 3. Bidirectional helper

```r
both <- selcorr_bidir(x, y, selection = selection_spec("de_anchor"), method = "replay")
both$forward   # x anchored
both$reverse   # y anchored
```

Two separate estimands. Null centres, dispersions and p values are never
averaged or combined.

## 4. Subordinate functions

### `diagnose(fit, type = "recovery", fractions, effect_sizes, n_outer, n_perm, seed, checkpoint, parallel, alpha)`

Injects direction-aligned shared effects into the replayed arm of a
sample-level fit and reports the rejection rate of the replay test per grid
cell. Returns `fraction`, `effect_size`, `rejection_rate`, `mc_se`,
`median_observed_r`, `median_null_mean`, `median_p`.

### `validate(fit, x_new, y_new, n_perm, seed, parallel)`

Evaluates the selected gene set in data that took no part in selection. The set
is externally fixed and the evaluation reference does not replay the original
selection. Returns `n_eval_common`, `n_eval_selected`, `evaluation_r`,
`null_mean`, `null_sd`, `critical_value`, `percentile`, `p_value`. The attached
note records that this p value is not combined with the replay p value.

## 5. Estimator adapters

```r
estimator_limma(expr, design_full, coef_pattern = "^group")
estimator_deseq2(counts, colData, design = ~group, contrast = c("group", "case", "control"))
```

`limma` uses `trend = TRUE` and `robust = TRUE`, as in the archived workflows.
`DESeq2` is provided for completeness; substituting an estimator is an analysis
choice and must be reported as such.

## 6. Behavioural contract

* A draw whose selected set falls below `min_selected` yields `NA` and is
  counted as failed, never silently dropped from the denominator bookkeeping.
* The fixed-set comparator is always flagged `diagnostic_only = TRUE`.
* Draw sequences depend only on `seed`, `B` and the arm sizes.
* Checkpoints let an interrupted run resume; completed draws are reused
  unchanged.
* No function in the package asserts that a reference is unbiased, universally
  valid, or that a non-significant result proves absence of biological sharing.

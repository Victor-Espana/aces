<table style="border:0;width:100%;border-collapse:collapse;margin:0;">
  <tr>
    <td style="width:200px;vertical-align:middle;padding:0;border:0;">
      <img src="man/figures/aces_logo.png" alt="aces" width="190" style="display:block;margin:0;">
    </td>
    <td style="vertical-align:middle;padding:0 0 0 10px;border:0;">
      <h1 style="margin:0;font-weight:700;font-size:3.4rem;line-height:1;color:#4b5563;">
        Adaptive Constrained Enveloping Splines
      </h1>
    </td>
  </tr>
</table>
<hr style="margin:8px 0;border:0;border-top:1px solid #e5e7eb;">

# Overview

`aces` estimates maximum attainable outputs with adaptive regression splines,
uses them to construct production technologies, and measures technical
efficiency. It provides three related estimators:

- **ACES** is the standard estimator.
- **RF-ACES** fits an ensemble of randomized ACES learners.
- **Quick ACES** reduces the candidate search for faster fitting.

The package also computes radial, directional, Russell, and weighted additive
efficiency measures and their efficient targets.

## Choosing an estimator

| Estimator | Recommended use | Main trade-off |
|:--|:--|:--|
| ACES | Maximum-output and technology estimation | Fits and prunes a full set of candidate basis functions |
| RF-ACES | Greater stability through an ensemble of randomized models | Usually requires more computation |
| Quick ACES | Exploratory work or larger candidate grids | Gains speed by considering fewer candidates |

## How ACES builds the technology

ACES separates output estimation from technology construction:

1. **Estimate maximum outputs.** The forward step builds a flexible spline model
   shared by all outputs, and the backward step removes terms to control
   complexity. If `shape$mono` and `shape$conc` are both `TRUE`, ACES estimates a
   non-decreasing and concave production function for each output. If either is
   `FALSE`, the spline model instead pushes observed outputs toward their
   estimated maximum production capacity.
2. **Refine multi-output predictions.** Separate output models do not necessarily
   produce a jointly feasible output vector. ACES compares joint and
   output-specific radial scores and retains observed components when a marginal
   prediction may overstate joint production capacity. With one output, all
   reference outputs come from the fitted model.
3. **Construct the technology.** The original inputs and refined output vectors
   form the reference sample. Efficiency scores and targets apply a DEA-type
   envelopment to these points under the selected returns-to-scale assumption.

A production technology is therefore constructed for every setting of `shape`.
The argument controls the interpretation of the spline stage, not whether the
final technology exists.

# Installation

```r
# From CRAN
install.packages("aces")

# Development version from GitHub
# install.packages("devtools")
devtools::install_github("Victor-Espana/aces")
```

# Basic workflow

```r
library(aces)

set.seed(314)
production_data <- cobb_douglas_XnY1(N = 50, nX = 3)

fit <- aces(
  data = production_data,
  x = 1:3,
  y = 4,
  mul_BF = list(
    max_degree = 2,
    inter_cost = 0.05
  )
)

# Points that define the estimated technology
technology <- get_technology(fit, method = "aces")

# Output-oriented efficiency scores
scores <- get_scores(
  eval_data = production_data,
  x = 1:3,
  y = 4,
  object = fit,
  measure = "rad_out"
)

# Efficient input and output targets
targets <- get_targets(
  eval_data = production_data,
  x = 1:3,
  y = 4,
  object = fit,
  measure = "rad_out"
)
```

The fitted object stores the forward model, the selected pruned model, optional
cubic and quintic smooth versions, and the technology associated with each one.
Use the same `method` in `get_technology()`, `get_scores()`, and `get_targets()`
when the results must refer to the same estimated technology.

# Main controls

The following arguments have the greatest effect on model complexity:

- `mul_BF$max_degree` sets the highest interaction degree. A value of 1 gives an
  additive model; larger values allow interactions between inputs. Prepared
  interactions use statistical `:` notation, for example `capital:labour`.
- `mul_BF$inter_cost` sets the extra improvement required before an interaction
  is selected. Larger values favor simpler models.
- `max_terms` limits the number of terms created during the forward step. Fitting
  may stop earlier when no candidate meets `err_red`.
- `err_red` is the minimum relative error reduction required to add another pair
  of first-degree basis functions. Larger values stop the search sooner.
- `shape$mono` and `shape$conc` control the spline estimation stage. When both
  are `TRUE`, the output-specific fits are production functions. If either is
  `FALSE`, they are intermediate predictors of maximum attainable output. ACES
  still constructs the final production technology from the refined output
  vectors.
- `kn_grid` supplies candidate knots. Leave it unspecified for the standard
  observed-value grid, or provide a list to control candidates by input.

# Efficiency measures

Choose a measure according to the adjustment that is meaningful for the
application:

| Measure | Adjustment represented |
|:--|:--|
| `rad_out` | Expands all outputs proportionally with inputs fixed |
| `rad_inp` | Contracts all inputs proportionally with outputs fixed |
| `ddf` | Contracts inputs and expands outputs along a chosen direction |
| `rsl_out`, `rsl_inp` | Allows a separate proportional adjustment for each output or input |
| `wam_*` | Uses input and output slacks with the selected weighting rule |

Set `returns = "variable"` for a convex technology with a convexity constraint,
or `returns = "constant"` for a conical constant-returns technology.

# RF-ACES

RF-ACES fits ACES learners on randomized samples and input subsets, then
aggregates their predictions and technologies. This can make the estimated
technology less sensitive to a single sample or model, but fitting the ensemble is
more expensive. `learners`, `bag_size`, and `max_feats` control its size and
randomization. The same interpretation of `shape` applies within every learner.

```r
rf_fit <- rf_aces(
  data = production_data,
  x = 1:3,
  y = 4,
  learners = 20,
  max_feats = 1
)

predicted_outputs <- predict(
  rf_fit,
  newdata = production_data,
  x = 1:3
)

# Empirical 95% intervals across RF-ACES learners
intervals <- rf_aces_intervals(
  rf_fit,
  newdata = production_data,
  x = 1:3,
  y = 4,
  level = 0.95,
  type = "both"
)

intervals$output
intervals$score

# Approximate OOB-calibrated prediction intervals for observed outputs
prediction_intervals <- rf_aces_intervals(
  rf_fit,
  newdata = production_data,
  x = 1:3,
  type = "output",
  level = 0.95,
  calibration = "oob"
)
```

These limits are empirical quantiles across the fitted learners. They describe
the forest's sensitivity to bootstrap sampling and randomized input selection;
they are not calibrated confidence or prediction intervals. Setting
`calibration = "oob"` instead uses out-of-bag residuals to construct approximate
prediction intervals for new observed outputs. It does not estimate coverage for
the latent production frontier or for true efficiency scores.

# Variable importance and relevant inputs

The package uses several related, but different, notions of input importance:

| Mechanism | What it measures | What it is used for |
|:--|:--|:--|
| Quick ACES correlation filter | Rank association between each prepared input and the output | Excludes weak variables from the forward candidate search |
| `aces_varimp(..., importance = "forward")` | Best training-error reduction attainable with each evaluated prepared input in every forward iteration | Produces an ACES ranking, or the learner-averaged ranking for RF-ACES |
| `aces_varimp(..., importance = "permutation")` | Loss increase after permuting one original input and rebuilding its interactions | Measures sensitivity of the fixed ACES or RF-ACES fit |
| `aces_varimp(..., importance = "loco")` | Loss increase after refitting without one original input and its interactions | Measures the input's contribution to the complete ACES or RF-ACES procedure |
| Selected basis functions | Whether an original input appears directly or through an interaction | Defines the reduced input set used when `relevant = TRUE` |

These quantities are not interchangeable. In particular, `relevant = TRUE` in
`get_scores()` and `get_targets()` does not apply a threshold to
`aces_varimp()`. It retains inputs represented in the basis functions selected
by the requested ACES model. For RF-ACES, an input is retained if it appears in
at least one learner. The reduced set changes the input dimensions of the
efficiency problem; it does not refit the spline model.

```r
# Structural selection based on the fitted basis functions
reduced_scores <- get_scores(
  eval_data = production_data,
  x = 1:3,
  y = 4,
  relevant = TRUE,
  object = fit
)

# OOB predictive ranking for RF-ACES; not used automatically by reduced_scores
importance <- aces_varimp(
  rf_fit,
  importance = "permutation",
  control = list(method = "rf_aces", eval_data = "oob")
)
```

# Quick ACES

Set `quick_aces = TRUE` to activate three reductions before and during the
forward search: rank-correlation filtering, DEA-guided knot neighborhoods around
efficient DMUs, and adaptive candidate allocation based on recent error
reductions. The correlation filter retains an input when its Spearman or Kendall
correlation reaches the corresponding threshold for at least one output. The
model retains the same output structure as standard ACES, but the smaller search
can select a different model. Its correlation and error-reduction values are
search diagnostics rather than test-set or causal importance measures. The
error-reduction history and its cumulative ranking are returned by
`aces_varimp()`.

```r
quick_fit <- aces(
  data = production_data,
  x = 1:3,
  y = 4,
  quick_aces = TRUE
)

# Cumulative reduction ranking and iteration-by-iteration history
aces_varimp(quick_fit)
aces_varimp(quick_fit, control = list(type = "iterations"))

# Fixed-model permutation and refitted leave-one-covariate-out importance
aces_varimp(
  quick_fit,
  importance = "permutation",
  control = list(repeats = 20, seed = 314)
)
aces_varimp(quick_fit, importance = "loco")
```

# Citing

Please cite the method that you use:

- España, V. J., Aparicio, J., Barber, X., & Esteve, M. (2024). *Estimating
  production functions through additive models based on regression splines*.
  European Journal of Operational Research, 312(2), 684–699.
- España, V. J., Aparicio, J., & Barber, X. (2025). *Estimating production
  technologies using multi-output adaptive constrained enveloping splines*.
  Computers & Operations Research, 107242.
- España, V. J., Aparicio, J., & Barber, X. (2025). *An adaptation of Random
  Forest to estimate convex non-parametric production technologies: an empirical
  illustration of efficiency measurement in education*. International
  Transactions in Operational Research, 32(5), 2523–2546.

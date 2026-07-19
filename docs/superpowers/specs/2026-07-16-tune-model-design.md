# Design: simple tuning API for ACES and RF-ACES

## Public API

```r
partition <- data_partition(
  data,
  x = 1:6,
  y = 7,
  validation = "kfold",
  resampling_args = list(folds = 5),
  delta = 1.5,
  seed = 42
)

aces_fit <- tune_model_aces(
  partition,
  mul_BF = list(
    max_degree = c(1, 2),
    inter_cost = c(0.01, 0.05)
  )
)

rf_fit <- tune_model_rf_aces(
  partition,
  mul_BF = list(max_degree = c(1, 2), inter_cost = 0.05),
  learners = c(50, 100),
  max_feats = c(1, 2)
)
```

`data_partition()` receives input and output column indexes. It owns the data,
resampling indexes, `delta`, and seed, so the same partition can be reused by
both tuning functions.

The model arguments in `tune_model_aces()` and `tune_model_rf_aces()` have the
same names, defaults, and nested form as in `aces()` and `rf_aces()`. Vectors
inside `mul_BF`, `shape`, and `early_stopping` are expanded internally.

## Internal flow

1. Build the Cartesian product of model parameters with `expand.grid()`.
2. Compute DEA validation weights for each stored split using the partition's
   `delta`.
3. Fit and evaluate every configuration on every split.
4. Select the configuration with the smallest mean validation error.
5. Fit that configuration once more using all observations.

The implementation uses one shared tuning loop and small base R helpers.

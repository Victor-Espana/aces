# ==============================================
# Test script: Random Forest ACES functions
# ==============================================

devtools::load_all()

set.seed(42)

# --- 1. Simulate data ---
data <- cobb_douglas_XnY1(N = 100, nX = 3)
x <- 1:3
y <- 4

cat("=== Data generated: 100 obs, 3 inputs, 1 output ===\n")
cat("Head:\n")
print(head(data[, c(x, y)]))

# --- 2. Train RF-ACES ---
cat("\n=== Training RF-ACES (10 trees) ===\n")
model <- rf_aces(
  data = data,
  x = x,
  y = y,
  mul_BF = list(max_degree = 2, inter_cost = 0.05),
  learners = 100,
  metric = "mae",
  early_stopping = list("ma_window" = 10, "tolerance" = 5)
)

cat("\nClass:", class(model), "\n")
cat("Trees trained:", length(model[["forest"]]), "\n")
cat("Metric:", model[["control"]][["metric"]], "\n")

# --- 3. Check sample_bag stored ---
cat("\n=== Checking sample_bag storage ===\n")
has_bags <- all(sapply(model[["forest"]], function(t) !is.null(t[["sample_bag"]])))
cat("All trees have sample_bag:", has_bags, "\n")
cat("Bag size (tree 1):", length(model[["forest"]][[1]][["sample_bag"]]), "\n")

# --- 4. Predict (S3 method) ---
cat("\n=== predict.rf_aces ===\n")
preds <- predict(model, newdata = data, x = x)
cat("Predictions shape:", nrow(preds), "x", ncol(preds), "\n")
cat("First 5 predictions:\n")
print(head(preds, 5))

# --- 5. Variable importance ---
cat("\n=== rf_aces_varimp ===\n")
varimp <- rf_aces_varimp(
  data = data,
  x = x,
  y = y,
  object = model,
  repeats = 2
)
cat("Variable importance:\n")
print(varimp)

# --- 6. get_scores (unified) with rf_aces object ---
cat("\n=== get_scores (rad_out) ===\n")
scores_rad <- get_scores(
  eval_data = data,
  x = x,
  y = y,
  object = model,
  measure = "rad_out"
)
cat("Scores shape:", nrow(scores_rad), "x", ncol(scores_rad), "\n")
cat("Summary:\n")
print(summary(scores_rad))

# --- 7. get_scores with rf_aces_rad_out ---
cat("\n=== get_scores (rf_aces_rad_out) ===\n")
scores_rf <- get_scores(
  eval_data = data,
  x = x,
  y = y,
  object = model,
  measure = "rf_aces_rad_out"
)
cat("RF-ACES rad_out scores:\n")
print(summary(scores_rf))

# --- 8. get_targets (unified) ---
cat("\n=== get_targets (rad_out) ===\n")
targets <- get_targets(
  eval_data = data,
  x = x,
  y = y,
  object = model,
  measure = "rad_out"
)
cat("Targets shape:", nrow(targets), "x", ncol(targets), "\n")
cat("First 5 targets:\n")
print(head(targets, 5))

# --- 9. get_targets with rf_aces_rad_out ---
cat("\n=== get_targets (rf_aces_rad_out) ===\n")
targets_rf <- get_targets(
  eval_data = data,
  x = x,
  y = y,
  object = model,
  measure = "rf_aces_rad_out"
)
cat("RF-ACES targets:\n")
print(head(targets_rf, 5))

# --- 10. get_scores with relevant = TRUE ---
cat("\n=== get_scores (relevant = TRUE) ===\n")
scores_rel <- get_scores(
  eval_data = data,
  x = x,
  y = y,
  relevant = TRUE,
  object = model,
  measure = "rad_out"
)
cat("Relevant scores:\n")
print(summary(scores_rel))

cat("\n=== ALL TESTS PASSED ===\n")

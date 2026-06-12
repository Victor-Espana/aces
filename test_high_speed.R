# ==================================================================== #
#  TEST DEL MODO HIGH-SPEED — paquete `aces`                            #
# ==================================================================== #
#
# UNA SOLA EJECUCION, SIN CAMBIAR DE RAMA:
#
#   source("test_high_speed.R")
#
# Compara, en la misma sesion, aces(high_speed = FALSE) [implementacion
# original, intacta] contra aces(high_speed = TRUE) [implementaciones
# *_high_speed], y verifica:
#
#   PARTE A — Equivalencia matematica (criterio de PASS/FAIL):
#     A1. Scores FDH: rad_out_fdh (enumeracion) vs MILP legacy.
#     A2. LP de coeficientes: en CADA llamada real durante varios ajustes
#         se ejecutan las dos implementaciones y se comparan valor
#         objetivo y predicciones. (Los coeficientes pueden diferir si el
#         LP tiene optimos multiples: vertices distintos e igualmente
#         optimos.)
#     A3. Lo mismo para el LP de suavizado (cubico/quintico).
#
#   PARTE B — Comparacion end-to-end y tiempos (informativa):
#     Ajusta cada escenario con ambos modos y compara nudos,
#     coeficientes, predicciones y scores, ademas de medir tiempos.
#     Diferencias aqui NO implican error: si A1-A3 pasan, ambos modelos
#     son igualmente optimos paso a paso; los caminos pueden divergir en
#     empates exactos entre nudos candidatos.
#
# Duracion: la parte legacy de los escenarios grandes puede tardar
# varios minutos.
# ==================================================================== #

suppressMessages({
  devtools::load_all(".", quiet = TRUE)
})

TOL <- 1e-6

run_quiet <- function(expr) {
  out <- capture.output(res <- force(expr))
  res
}

cat("\n==============================================================\n")
cat("  TEST HIGH-SPEED — una sola ejecucion\n")
cat("==============================================================\n")

fallos_A <- 0

# ==================================================================== #
#  PARTE A1 — FDH: enumeracion vs MILP legacy                           #
# ==================================================================== #

cat("\n--- A1: scores FDH (rad_out_fdh vs MILP legacy) ---\n")

a1 <- list()

for (cfg in list(c(N = 60, nX = 1, nY = 1),
                 c(N = 80, nX = 3, nY = 1),
                 c(N = 60, nX = 2, nY = 2),
                 c(N = 70, nX = 3, nY = 3))) {
  set.seed(1000 + cfg[["N"]])
  xmat <- matrix(runif(cfg[["N"]] * cfg[["nX"]], 1, 10), ncol = cfg[["nX"]])
  ymat <- matrix(runif(cfg[["N"]] * cfg[["nY"]], 1, 10), ncol = cfg[["nY"]])

  # legacy: MILP con binarias (via opcion en FALSE)
  options(aces.high_speed = FALSE)
  s_leg <- aces:::rad_out(
    tech_xmat = xmat, tech_ymat = ymat,
    eval_xmat = xmat, eval_ymat = ymat,
    convexity = FALSE, returns = "variable", type = "objective"
  )
  options(aces.high_speed = TRUE)

  # high-speed: enumeracion
  s_hs <- aces:::rad_out_fdh(
    tech_xmat = xmat, tech_ymat = ymat,
    eval_xmat = xmat, eval_ymat = ymat
  )

  dmax <- max(abs(s_leg - s_hs))
  estado <- if (dmax <= TOL) "PASS" else { fallos_A <- fallos_A + 1; "FAIL" }
  cat(sprintf("  [%s] N=%d nX=%d nY=%d  | dif. max = %.2e\n",
              estado, cfg[["N"]], cfg[["nX"]], cfg[["nY"]], dmax))
  a1[[paste0("N", cfg[["N"]], "_x", cfg[["nX"]], "_y", cfg[["nY"]])]] <- dmax
}

# ==================================================================== #
#  PARTE A2 / A3 — LPs: ambas implementaciones en cada llamada real     #
# ==================================================================== #

cat("\n--- A2: LP de coeficientes (lineal) | A3: LP de suavizado ---\n")

audit <- new.env()
audit$n_lin <- 0L; audit$obj_lin <- 0; audit$yhat_lin <- 0; audit$coef_lin <- 0
audit$n_smo <- 0L; audit$obj_smo <- 0; audit$yhat_smo <- 0; audit$coef_smo <- 0
audit$fallback <- 0L

orig_ec  <- get("estimate_coefficients", envir = asNamespace("aces"))
orig_ecs <- get("estimate_coefficients_smoothed", envir = asNamespace("aces"))

audit_ec <- function(B, y_obs, dea_scores, fdh_scores, it_list, Bp_list, shape) {
  c_hs <- aces:::estimate_coefficients_high_speed(
    B = B, y_obs = y_obs, dea_scores = dea_scores, fdh_scores = fdh_scores,
    it_list = it_list, Bp_list = Bp_list, shape = shape
  )
  c_lg <- aces:::estimate_coefficients_envelopment(
    B = B, y_obs = y_obs, dea_scores = dea_scores, fdh_scores = fdh_scores,
    it_list = it_list, Bp_list = Bp_list, shape = shape
  )
  deaw <- as.matrix(1 / dea_scores)
  for (out in 1:ncol(y_obs)) {
    env_ind <- fdh_scores[, out] < 1 + 1e-5
    w <- deaw[env_ind, out]
    B_env <- B[env_ind, , drop = FALSE]
    obj_h <- sum(w * (B_env %*% c_hs[, out] - y_obs[env_ind, out]))
    obj_l <- sum(w * (B_env %*% c_lg[, out] - y_obs[env_ind, out]))
    audit$obj_lin  <- max(audit$obj_lin, abs(obj_h - obj_l) / max(1, abs(obj_l)))
    audit$yhat_lin <- max(audit$yhat_lin, max(abs(B %*% (c_hs[, out] - c_lg[, out]))))
    audit$coef_lin <- max(audit$coef_lin, max(abs(c_hs[, out] - c_lg[, out])))
  }
  audit$n_lin <- audit$n_lin + 1L
  c_hs
}

audit_ecs <- function(B, y_obs, dea_scores, n_pair, n_lsub, shape) {
  c_hs <- aces:::estim_coefs_smooth_high_speed(
    B = B, y_obs = y_obs, dea_scores = dea_scores,
    n_pair = n_pair, n_lsub = n_lsub, shape = shape
  )
  c_lg <- aces:::estim_coefs_smooth_envelopment(
    B = B, y_obs = y_obs, dea_scores = dea_scores,
    n_pair = n_pair, n_lsub = n_lsub, shape = shape
  )
  deaw <- as.matrix(1 / dea_scores)
  for (out in 1:ncol(y_obs)) {
    w <- deaw[, out]
    obj_h <- sum(w * (B %*% c_hs[, out] - y_obs[, out]))
    obj_l <- sum(w * (B %*% c_lg[, out] - y_obs[, out]))
    audit$obj_smo  <- max(audit$obj_smo, abs(obj_h - obj_l) / max(1, abs(obj_l)))
    audit$yhat_smo <- max(audit$yhat_smo, max(abs(B %*% (c_hs[, out] - c_lg[, out]))))
    audit$coef_smo <- max(audit$coef_smo, max(abs(c_hs[, out] - c_lg[, out])))
  }
  audit$n_smo <- audit$n_smo + 1L
  c_hs
}

assignInNamespace("estimate_coefficients", audit_ec, ns = "aces")
assignInNamespace("estimate_coefficients_smoothed", audit_ecs, ns = "aces")

set.seed(123)
d1 <- cobb_douglas_XnY1(60, 1)
set.seed(123); invisible(run_quiet(aces(data = d1, x = 1, y = 2)))

set.seed(456)
d2 <- cobb_douglas_XnY1(60, 3)
set.seed(456); invisible(run_quiet(aces(data = d2, x = 1:3, y = 4)))

set.seed(789)
d3 <- cobb_douglas_X3Y3(50, border = 0.2, noise = 0, returns = "CRS")
set.seed(789); invisible(run_quiet(aces(data = d3, x = 1:3, y = 4:6)))

assignInNamespace("estimate_coefficients", orig_ec, ns = "aces")
assignInNamespace("estimate_coefficients_smoothed", orig_ecs, ns = "aces")

for (bloque in list(
  list(et = "A2 lineal ", n = audit$n_lin, o = audit$obj_lin,
       y = audit$yhat_lin, co = audit$coef_lin),
  list(et = "A3 suavizado", n = audit$n_smo, o = audit$obj_smo,
       y = audit$yhat_smo, co = audit$coef_smo)
)) {
  ok <- bloque$o <= 1e-7 && bloque$y <= TOL
  if (!ok) fallos_A <- fallos_A + 1
  cat(sprintf("  [%s] %s: %d LPs | dif.rel.objetivo = %.2e | dif.predicciones = %.2e | dif.coefs = %.2e\n",
              if (ok) "PASS" else "FAIL", bloque$et, bloque$n, bloque$o, bloque$y, bloque$co))
}
cat("  (dif.coefs > 0 con objetivo y predicciones identicos = optimos multiples; no es un error)\n")

# ==================================================================== #
#  PARTE B — end-to-end: high_speed = FALSE vs TRUE                     #
# ==================================================================== #

cat("\n--- B: escenarios end-to-end (legacy vs high-speed) ---\n")

escenarios <- list(
  esc1_1x_N100       = list(gen = function() cobb_douglas_XnY1(100, 1),
                            x = 1,   y = 2,   quick = FALSE, seed = 11),
  esc2_3x_N150       = list(gen = function() cobb_douglas_XnY1(150, 3),
                            x = 1:3, y = 4,   quick = FALSE, seed = 22),
  esc3_3x3y_N100     = list(gen = function() cobb_douglas_X3Y3(100, border = 0.1,
                                                               noise = 0, returns = "CRS"),
                            x = 1:3, y = 4:6, quick = FALSE, seed = 33),
  esc4_1x_N500       = list(gen = function() cobb_douglas_XnY1(500, 1),
                            x = 1,   y = 2,   quick = FALSE, seed = 44),
  esc5_3x_N300       = list(gen = function() cobb_douglas_XnY1(300, 3),
                            x = 1:3, y = 4,   quick = FALSE, seed = 55),
  esc6_3x_N300_quick = list(gen = function() cobb_douglas_XnY1(300, 3),
                            x = 1:3, y = 4,   quick = TRUE,  seed = 66)
)

extraer <- function(m) {
  out <- list()
  for (sub in names(m[["methods"]])) {
    sm <- m[["methods"]][[sub]]
    if (is.null(sm)) next
    el <- list()
    if (!is.null(sm[["knots"]])) el[["knots"]] <- sm[["knots"]]
    if (!is.null(sm[["coefs"]])) el[["coefs"]] <- as.matrix(sm[["coefs"]])
    if (!is.null(sm[["Bmatx"]]) && !is.null(sm[["coefs"]])) {
      el[["yhat"]] <- as.matrix(sm[["Bmatx"]]) %*% as.matrix(sm[["coefs"]])
    }
    out[[sub]] <- el
  }
  out
}

max_diff <- function(a, b) {
  if (is.null(a) && is.null(b)) return(0)
  if (is.null(a) || is.null(b)) return(Inf)
  a <- as.matrix(a); b <- as.matrix(b)
  if (!all(dim(a) == dim(b))) return(Inf)
  suppressWarnings(max(abs(a - b)))
}

knots_iguales <- function(ka, kb) {
  if (is.null(ka) && is.null(kb)) return(TRUE)
  if (is.null(ka) || is.null(kb)) return(FALSE)
  if (!all(dim(ka) == dim(kb))) return(FALSE)
  identical(ka$xi, kb$xi) && identical(ka$status, kb$status) &&
    max(abs(ka$t - kb$t)) <= TOL
}

resultados <- list()

cat(sprintf("\n  %-20s %10s %10s %8s | %-10s %-12s %-12s\n",
            "escenario", "legacy(s)", "high(s)", "speedup",
            "nudos", "dif.yhat", "dif.scores"))

for (nombre in names(escenarios)) {
  esc <- escenarios[[nombre]]
  set.seed(esc$seed)
  df <- esc$gen()

  res <- list()

  for (modo in c("legacy", "high")) {
    hs <- (modo == "high")
    set.seed(esc$seed)
    t_fit <- system.time(
      m <- run_quiet(aces(data = df, x = esc$x, y = esc$y,
                          quick_aces = esc$quick, high_speed = hs))
    )[["elapsed"]]
    s <- run_quiet(get_scores(eval_data = df, x = esc$x, y = esc$y,
                              object = m, method = "aces",
                              measure = "rad_out"))
    res[[modo]] <- list(t = t_fit, scores = as.matrix(s), sub = extraer(m))
  }

  # comparacion
  kn_ok <- TRUE; d_yhat <- 0; d_coef <- 0
  subs <- union(names(res$legacy$sub), names(res$high$sub))
  for (sub in subs) {
    a <- res$legacy$sub[[sub]]; b <- res$high$sub[[sub]]
    kn_ok  <- kn_ok && knots_iguales(a$knots, b$knots)
    d_yhat <- max(d_yhat, max_diff(a$yhat, b$yhat))
    d_coef <- max(d_coef, max_diff(a$coefs, b$coefs))
  }
  d_sc <- max_diff(res$legacy$scores, res$high$scores)
  sp <- res$legacy$t / max(res$high$t, 1e-9)

  cat(sprintf("  %-20s %10.2f %10.2f %7.1fx | %-10s %-12s %-12s\n",
              nombre, res$legacy$t, res$high$t, sp,
              if (kn_ok) "iguales" else "DISTINTOS",
              sprintf("%.2e", d_yhat), sprintf("%.2e", d_sc)))

  resultados[[nombre]] <- list(
    t_legacy = res$legacy$t, t_high = res$high$t, speedup = sp,
    kn_ok = kn_ok, d_yhat = d_yhat, d_coef = d_coef, d_scores = d_sc
  )
}

# ==================================================================== #
#  VEREDICTO                                                           #
# ==================================================================== #

sps <- sapply(resultados, `[[`, "speedup")
identicos <- all(sapply(resultados, function(r) r$kn_ok && r$d_yhat <= TOL && r$d_scores <= TOL))

cat("\n==============================================================\n")
cat("  VEREDICTO\n")
cat("==============================================================\n")
cat(sprintf("  Equivalencia matematica (A1-A3): %s\n",
            if (fallos_A == 0) "PASS" else paste("FAIL en", fallos_A, "bloque(s)")))
cat(sprintf("  Speedup medio (geometrico): %.1fx | min: %.1fx | max: %.1fx\n",
            exp(mean(log(sps))), min(sps), max(sps)))
if (identicos) {
  cat("  Modelos finales: IDENTICOS en todos los escenarios.\n")
} else {
  cat("  Modelos finales: difieren en algun escenario. Si A1-A3 son PASS,\n")
  cat("  la causa son optimos multiples del LP: en empates exactos entre\n")
  cat("  nudos candidatos cada modo puede elegir un camino distinto e\n")
  cat("  igualmente optimo. Ver informe_optimizacion.md, seccion 6.\n")
}

saveRDS(
  list(meta = list(fecha = Sys.time(), R = R.version.string,
                   a1 = a1,
                   a2 = list(n = audit$n_lin, obj = audit$obj_lin,
                             yhat = audit$yhat_lin, coef = audit$coef_lin),
                   a3 = list(n = audit$n_smo, obj = audit$obj_smo,
                             yhat = audit$yhat_smo, coef = audit$coef_smo)),
       escenarios = resultados),
  "test_high_speed_resultados.rds"
)
cat("\nResultados guardados en: test_high_speed_resultados.rds\n")

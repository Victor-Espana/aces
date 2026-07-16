# ===========================================================================
#  Test script: verificar que ambas versiones de aces funcionan correctamente
# ===========================================================================
#
#  Este script testea ambos paquetes (aces y aces_slow) si están instalados.
#  Genera un único dataset y ejecuta el pipeline completo en cada uno.
#
# --- Instrucciones de uso -------------------------------------------------
#
#  PASO 1: Instalar ambos paquetes (solo una vez)
#
#    # Primero, instalar aces_slow (versión legacy):
#    #   Abrir una terminal y hacer:
#    #     cd "ruta/al/paquete"
#    #     git checkout legacy
#    #   En R:
#    devtools::install()
#
#    # Después, instalar aces (versión optimizada):
#    #   En la terminal:
#    #     git checkout v2
#    #   En R:
#    devtools::install()
#
#    # A partir de aquí, ambos paquetes están instalados y se pueden
#    # usar de forma independiente con library(aces) o library(aces_slow).
#
#  PASO 2: Ejecutar este script
#
#    source("test/test_ambas_versiones.R")
#
# --------------------------------------------------------------------------

# ===========================================================================
#  Función auxiliar: ejecuta el pipeline completo sobre un dataset dado
# ===========================================================================
run_test <- function(pkg_name) {

  cat("\n")
  cat("=====================================================\n")
  cat("  Testeando paquete:", pkg_name, "\n")
  cat("  Versión:", as.character(packageVersion(pkg_name)), "\n")
  cat("=====================================================\n\n")

  # --- 1. Generar datos ---------------------------------------------------

  set.seed(42)

  # Ambos paquetes exportan las mismas funciones con los mismos nombres
  data <- cobb_douglas_XnY1(N = 50, nX = 3)

  x <- 1:3    # x1, x2, x3
  y <- 4      # y

  cat("Datos: Cobb-Douglas, 3 inputs, 1 output, N =", nrow(data), "\n\n")

  # --- 2. Ajustar modelo aces ---------------------------------------------

  cat("Ajustando modelo aces...\n\n")

  t_start <- Sys.time()

  modelo <- aces(
    data = data,
    x    = x,
    y    = y,
    scale_data = TRUE,
    max_terms  = 20,
    err_red    = 0.01
  )

  t_end <- Sys.time()
  elapsed <- round(difftime(t_end, t_start, units = "secs"), 2)

  cat("\nModelo ajustado en", elapsed, "segundos.\n\n")

  # --- 3. Imprimir el modelo -----------------------------------------------

  cat("--- Modelo ACES (backward) ---\n")
  print(modelo, method = "aces")
  cat("\n")

  # --- 4. Predicciones (targets) -------------------------------------------

  cat("Calculando targets (rad_out, VRS)...\n")

  targets <- get_targets(
    eval_data = data,
    x         = x,
    y         = y,
    object    = modelo,
    method    = "aces",
    measure   = "rad_out",
    returns   = "variable"
  )

  cat("  Primeras 5 filas de targets:\n")
  print(head(targets, 5))
  cat("\n")

  # --- 5. RMSE frente a la frontera teórica --------------------------------

  y_aces <- targets[, ncol(targets)]
  y_true <- data[, "yT"]
  rmse   <- sqrt(mean((y_aces - y_true)^2))

  cat("  RMSE (ACES vs frontera teórica):", round(rmse, 4), "\n")
  cat("  Tiempo de ajuste:", elapsed, "segundos\n")

  return(list(rmse = rmse, tiempo = as.numeric(elapsed)))
}

# ===========================================================================
#  Ejecutar tests
# ===========================================================================

resultados <- list()

# --- Test aces (versión optimizada) ----------------------------------------
if (requireNamespace("aces", quietly = TRUE)) {
  library(aces)
  resultados[["aces"]] <- run_test("aces")
  detach("package:aces", unload = TRUE)
} else {
  cat("\n[!] Paquete 'aces' no instalado. Saltando...\n")
  cat("    Para instalarlo: git checkout v2 && devtools::install()\n")
}

# --- Test aces_slow (versión legacy) ----------------------------------------
if (requireNamespace("aces_slow", quietly = TRUE)) {
  library(aces_slow)
  resultados[["aces_slow"]] <- run_test("aces_slow")
  detach("package:aces_slow", unload = TRUE)
} else {
  cat("\n[!] Paquete 'aces_slow' no instalado. Saltando...\n")
  cat("    Para instalarlo: git checkout legacy && devtools::install()\n")
}

# --- Comparación final -----------------------------------------------------

cat("\n")
cat("=====================================================\n")
cat("  RESUMEN COMPARATIVO\n")
cat("=====================================================\n\n")

if (length(resultados) == 2) {
  cat(sprintf("  %-12s  RMSE = %.4f   Tiempo = %.2f s\n",
              "aces",      resultados$aces$rmse,      resultados$aces$tiempo))
  cat(sprintf("  %-12s  RMSE = %.4f   Tiempo = %.2f s\n",
              "aces_slow", resultados$aces_slow$rmse, resultados$aces_slow$tiempo))

  cat("\n")
  speedup <- resultados$aces_slow$tiempo / max(resultados$aces$tiempo, 0.01)
  cat(sprintf("  Speedup (aces_slow / aces): %.1fx\n", speedup))

  rmse_diff <- abs(resultados$aces$rmse - resultados$aces_slow$rmse)
  cat(sprintf("  Diferencia absoluta RMSE:   %.6f\n", rmse_diff))

  if (rmse_diff < 0.5) {
    cat("\n  [OK] Ambos paquetes producen resultados comparables.\n")
  } else {
    cat("\n  [!] Diferencia notable en RMSE. Revisar resultados.\n")
  }
} else if (length(resultados) == 1) {
  pkg <- names(resultados)[1]
  cat(sprintf("  Solo se pudo testear '%s'. Instala ambos para comparar.\n", pkg))
} else {
  cat("  No se pudo testear ningún paquete.\n")
}

cat("\n")

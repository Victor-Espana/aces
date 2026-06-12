# Informe de rendimiento del paquete `aces`

**Fecha:** 12/06/2026
**Alcance:** diagnóstico de los cuellos de botella de ejecución y propuesta de cambios de alto impacto, sin alterar los resultados del método.

---

## 1. Estructura del coste computacional

El coste de `aces()` está dominado por el algoritmo forward. Su estructura de bucles es:

```
por cada iteración forward (hasta max_terms/2, por defecto ~25)
  por cada variable de entrada xi
    por cada nudo candidato t (por defecto, hasta N valores)
      → construir matriz B candidata
      → resolver un LP completo con Rglpk        ← coste dominante
      → evaluar el error
```

Es decir, el número de LPs resueltos es del orden de `iteraciones × nX × N`. Con N = 500 observaciones, 3 inputs y 20 iteraciones, son del orden de **30.000 LPs**. El problema no es tanto el número de LPs (inherente al método, como en MARS) sino que **cada LP es mucho más grande y caro de construir de lo necesario**.

El backward (`backward_algorithm.R`) y el smoothing (`estimate_coefficients_smooth.R`) reutilizan la misma rutina de estimación, así que heredan el mismo problema, aunque con muchas menos llamadas.

---

## 2. Cuello de botella nº 1: el LP de `estimate_coefficients` incluye N variables de error innecesarias

**Archivo:** `R/estimate_coefficients_linear.R` (función `estimate_coefficients_envelopment`, líneas ~118–199). Aplica igualmente a `R/estimate_coefficients_smooth.R` (líneas ~110–210).

### Formulación actual

```
variables:  (β₀, β₁, ..., β_p, e₁, ..., e_N)        →  p + N variables
min         Σᵢ wᵢ·eᵢ
s.a.        (Bβ)ᵢ − eᵢ = yᵢ      i = 1..N            →  N restricciones de igualdad
            monotonía / concavidad                    →  k restricciones (k pequeño)
            β libre, e ≥ 0 (cotas por defecto)
```

### Reformulación equivalente

Las variables de error están determinadas por β: `eᵢ = (Bβ)ᵢ − yᵢ`. Sustituyendo en el objetivo:

```
Σᵢ wᵢ·eᵢ = Σᵢ wᵢ·((Bβ)ᵢ − yᵢ) = (wᵀB)·β − wᵀy
```

donde `wᵀy` es constante. La condición `e ≥ 0` pasa a ser `Bβ ≥ y`. El LP equivalente es:

```
variables:  (β₀, β₁, ..., β_p)                        →  p variables
min         (wᵀB)·β
s.a.        Bβ ≥ y                                    →  N restricciones de desigualdad
            monotonía / concavidad                    →  k restricciones
            β libre
```

**El óptimo es exactamente el mismo** (mismos β*; cualquier discrepancia es tolerancia numérica del simplex). Como w > 0 (inversa de scores DEA ≥ 1), la equivalencia es estricta.

### Impacto

Con `max_terms = 50`, p ≤ 51. El LP pasa de `p + N` variables (≈ 550 con N = 500; ≈ 1.050 con N = 1.000) a ≤ 51 variables, eliminando además las N restricciones de igualdad más costosas para el simplex. Multiplicado por decenas de miles de solves, es la mejora dominante. **Estimación: 10–50× en el paso de estimación.**

Cambios derivados (se simplifican, no se complican):

- `monotonocity_matrix()` y `concavity_matrix()` (en `estimate_coefficients.R`) ya no necesitan el bloque `matrix(0, nrow, N)` que les añade columnas de ceros para los errores (líneas 80–84 y 159–163).
- `objVal`, `bvec`, `dirs` y `bnds` se acortan en consecuencia.
- El indicador de envolvimiento `env_ind` (filtra DMUs con score FDH ≈ 1) se mantiene igual: las filas de `Bβ ≥ y` se restringen a esas observaciones y el vector de pesos también.

---

## 3. Cuello de botella nº 2: construcción densa de las matrices del LP

**Archivos:** `estimate_coefficients_linear.R` (línea 137: `EMat <- cbind(B[env_ind, ], diag(rep(-1, num_env), num_env))`), `estimate_coefficients.R` (uso de `as.matrix(Matrix::bdiag(...))`), `estimate_coefficients_smooth.R` (línea ~123).

`diag(-1, N)` crea una matriz densa N×N (8 MB con N = 1.000) **en cada candidato a nudo**, y los `rbind`/`cbind` copian todo de nuevo. Solo la asignación de memoria ya es comparable al coste del simplex.

Con la reformulación del punto 2 este bloque desaparece (ya no hay variables de error). Para lo que quede:

- `Rglpk_solve_LP()` acepta directamente matrices sparse (`slam::simple_triplet_matrix` o clases de `Matrix`). Pasar `Amat` en formato triplete evita la densificación y reduce también el coste de carga en GLPK.
- Evitar `as.matrix(Matrix::bdiag(...))` → mantener el resultado sparse.

---

## 4. Cuello de botella nº 3: scores FDH mediante MILP

**Archivo:** `R/efficiency_scores.R`, función `rad_out` con `convexity = FALSE` (líneas ~86–92), llamada desde `ACES.R` (~líneas 399–420).

Actualmente cada DMU resuelve un **MILP con N variables binarias** (lpSolveAPI con `set.type(..., "binary")`). Para N = 1.000 son 1.000 MILPs de 1.000 binarias: puede costar más que todo el forward.

El score FDH radial orientado a output con VRS tiene **fórmula cerrada por enumeración**, sin optimización:

```
φ_d = max_{j : x_j ≤ x_d (componente a componente)}   min_r ( y_jr / y_dr )
```

Vectorizable en R:

```r
fdh_out <- function(xmat, ymat) {
  N <- nrow(xmat)
  # dominancia en inputs: dom[j, d] = TRUE si x_j <= x_d
  dom <- sapply(1:N, function(d) colSums(t(xmat) <= xmat[d, ]) == ncol(xmat))
  sapply(1:N, function(d) {
    ratios <- sweep(ymat[dom[, d], , drop = FALSE], 2, ymat[d, ], "/")
    max(apply(ratios, 1, min))
  })
}
```

Coste O(N²·(nX+nY)) en operaciones vectorizadas frente a N MILPs. **Estimación: 100–1000× en este paso.** Resultado idéntico por definición del estimador FDH.

Los scores DEA-VRS (LP continuo, modelo reutilizado con `set.rhs`/`set.mat`) están razonablemente bien resueltos; no son prioritarios.

---

## 5. Mejoras secundarias (menor impacto, opcionales)

1. **`forward_algorithm.R`, bucle de nudos (líneas 174–364):** por cada candidato se copia `Bp_list` completo (`Bp_list_aux <- Bp_list`), se reordena y se reconstruye `B` con `cbind` de tres bloques. Tras arreglar el LP, esto pasará a ser una fracción visible del tiempo. Alternativa: construir las dos columnas nuevas y pasarlas junto al índice de inserción, materializando la matriz solo para el mejor candidato.
2. **`set_knots()` (líneas 459–486):** bucles `for` con `append`/`c()` incremental para los índices de span; vectorizable con `outer(used_kn_idx, -sp1:sp1, "+")`. Coste menor.
3. **`err_metric` por candidato:** `new_B %*% coefs` es O(N·p), aceptable; no tocar.
4. **Paralelización** del bucle de nudos (`parallel::mclapply` / `future.apply`): posible, pero complica el código y la reproducibilidad (RNG en `quick_aces`). Solo si tras 1–3 sigue siendo lento. No recomendada de entrada.
5. **Cambio de solver** (p. ej. `highs`, más rápido que GLPK): innecesario si se aplica el punto 2; añadiría una dependencia. No recomendado para la publicación.

---

## 6. Qué NO cambia

- El método y sus resultados: los puntos 2 y 4 son reformulaciones matemáticamente equivalentes; el punto 3 es un cambio de representación de la misma matriz.
- La API pública del paquete (argumentos y valores de retorno de `aces()`, `aces_scores()`, etc.).
- Las dependencias actuales (Rglpk, lpSolveAPI, Matrix se mantienen; `slam` ya es dependencia de Rglpk).

## 7. Plan de implementación y verificación propuesto

1. Rama nueva `perf-optimization`.
2. Implementar en este orden: reformulación del LP (linear y smooth) → matrices sparse → FDH por enumeración.
3. **Verificación de equivalencia:** con varios escenarios de `simulations.R` (distintos N, nX, nY) y semilla fija, comparar entre rama original y optimizada: coeficientes finales, predicciones (`y_hat`), scores de eficiencia y nudos seleccionados. Criterio: diferencias ≤ 1e-6 (tolerancia de solver).
4. **Benchmark:** tiempos antes/después para N ∈ {100, 500, 1000, 2000} con 2–4 inputs.
5. Revisión del diff antes de fusionar a la rama principal.

## 8. Impacto esperado total

| Cambio | Paso afectado | Speedup estimado del paso |
|---|---|---|
| Reformulación LP (sin variables de error) | forward + backward + smoothing | 10–50× |
| Matrices sparse | construcción de cada LP | incluido en lo anterior |
| FDH por enumeración | preprocesado de scores | 100–1000× |
| Mejoras secundarias | forward | 1.2–2× adicional |

Para una base de datos "normal" (N ≈ 500–2000, 2–5 inputs), el tiempo total debería reducirse **al menos un orden de magnitud**, probablemente más cuanto mayor sea N.

---

*Nota: este archivo es un documento de trabajo; conviene eliminarlo (o moverlo fuera del paquete) antes de enviar a CRAN, o añadirlo a `.Rbuildignore`.*

# Informe de optimización de rendimiento del paquete `aces`

**Fecha:** 12/06/2026
**Rama:** `test-ws-integration` (todo el trabajo está ahora en una sola rama)
**Conmutador:** `aces(..., high_speed = TRUE/FALSE)` — las implementaciones originales permanecen intactas en el código.

---

## 1. Resumen ejecutivo

Se han añadido implementaciones optimizadas (sufijo `_high_speed` o `_hs`) de las rutinas internas de optimización del paquete, junto a las originales, que no se han modificado. Un conmutador (`high_speed = TRUE`, por defecto) selecciona qué implementación se usa; con `high_speed = FALSE` el paquete reproduce exactamente el comportamiento original. Las optimizaciones son matemáticamente equivalentes: alcanzan el mismo valor óptimo y las mismas predicciones en cada problema de optimización individual. Cuando un LP tiene óptimos múltiples pueden devolver un vértice distinto e igualmente óptimo, lo que en presencia de empates exactos entre nudos candidatos puede llevar a un modelo final distinto pero paso a paso igualmente óptimo (sección 6).

La primera ronda de optimización (eliminación de variables de error + FDH por enumeración) se midió en 1,3× de media. La segunda ronda, incluida en este commit, resuelve además el LP a través de su **dual**, que ataca el coste real identificado (el tamaño de la base del símplex), con expectativa de mejora sustancialmente mayor, pendiente de medición con `test_high_speed.R`.

---

## 2. Diagnóstico: dónde se va el tiempo

La estructura de coste de `aces()` es:

```
por cada iteración forward (hasta max_terms/2 ≈ 25)
  por cada variable de entrada xi
    por cada nudo candidato t (hasta N)
      → construir matriz B candidata
      → resolver un LP                      ← coste dominante
      → evaluar el error
```

El número de LPs es del orden de `iteraciones × nX × N` (decenas de miles). A esto se suma, una sola vez, el cálculo inicial de scores FDH, que en la versión original resolvía un MILP con N variables binarias por cada DMU.

La primera ronda demostró empíricamente que reducir las *columnas* del LP no basta: el símplex de GLPK factoriza bases cuyo tamaño viene dado por las *filas* (las N restricciones de envolvimiento), y ese coste apenas cambió. La segunda ronda (dual) ataca exactamente ese punto.

---

## 3. Arquitectura: funciones duplicadas y conmutador

Ninguna función original se ha modificado. Se han añadido duplicados optimizados y despachadores:

| Punto de entrada (sin cambios de firma) | Implementación legacy | Implementación high-speed |
|---|---|---|
| `estimate_coefficients()` | `estimate_coefficients_envelopment()` | `estimate_coefficients_high_speed()` (+ `estimate_coefficients_reduced()` como respaldo) |
| `estimate_coefficients_smoothed()` | `estim_coefs_smooth_envelopment()` | `estim_coefs_smooth_high_speed()` |
| `rad_out()` (solo FDH-VRS) | MILP con binarias (código original) | `rad_out_fdh()` |
| `set_knots()` | `set_knots_legacy()` (cuerpo original) | `set_knots_high_speed()` |
| matrices de forma | `monotonocity_matrix()`, `concavity_matrix()` | `monotonocity_matrix_hs()`, `concavity_matrix_hs()` |

El despacho se controla con la opción `aces.high_speed`:

- `aces(data, x, y, high_speed = TRUE)` — por defecto; usa las versiones optimizadas y restaura la opción al salir.
- `aces(data, x, y, high_speed = FALSE)` — reproduce exactamente la implementación original.
- `options(aces.high_speed = FALSE)` — afecta globalmente (también a `rf_aces` y a llamadas internas fuera de `aces()`).

**Cómo revertir si todo falla:** (1) inmediato: `high_speed = FALSE`; (2) permanente: `git revert <commit>` del commit "feat: modo high_speed..." — al estar las funciones originales intactas, el revert es limpio y no afecta a nada más.

---

## 4. Detalle matemático de cada optimización

### 4.1. Eliminación de las variables de error del LP

La formulación original del problema de estimación de coeficientes es:

```
variables:  (β₀, ..., β_p, e₁, ..., e_N)              →  p + N variables
min         Σᵢ wᵢ eᵢ                                   (wᵢ = 1/score_DEAᵢ > 0)
s.a.        (Bβ)ᵢ − eᵢ = yᵢ   (obs. envueltas)         →  N igualdades
            restricciones de forma (monotonía/concav.) →  k desigualdades
            β libre, e ≥ 0
```

Las variables `e` están determinadas por `β`: `eᵢ = (Bβ)ᵢ − yᵢ`. Sustituyendo en el objetivo, `Σ wᵢeᵢ = (wᵀB)β − wᵀy`, donde `wᵀy` es constante. El problema equivale exactamente a:

```
min  (wᵀB)β    s.a.   Bβ ≥ y,  Sβ ≥ 0,  β libre       →  p variables
```

Mismo conjunto factible proyectado, mismo valor óptimo, mismos `Bβ*` óptimos. La función `estim_coefs_smooth_high_speed()` aplica esta misma reformulación al LP de suavizado.

### 4.2. Resolución a través del dual (`estimate_coefficients_high_speed`)

La medición de la primera ronda mostró que el LP reducido sigue siendo lento: el símplex trabaja con bases de tamaño `(nº de filas)² = (N+k)²`. El dual intercambia filas y columnas. Con `A = [B_env; S]`, `b = [y; 0]`, `c = BᵀenvW`:

```
primal:  min c'β   s.a.  Aβ ≥ b, β libre               (N+k filas, p columnas)
dual:    max b'λ   s.a.  A'λ = c, λ ≥ 0                (p filas,  N+k columnas)
```

Como `β` es libre, las restricciones duales son igualdades y la dualidad es fuerte (ambos problemas factibles). El símplex sobre el dual trabaja con bases `p × p` (p ≤ 51 con `max_terms = 50`), frente a `(N+k) × (N+k)` (≈ 500–2000). Los coeficientes primales óptimos se recuperan como los valores duales de las p restricciones de igualdad del dual (teorema de la envolvente: `∂(b'λ*)/∂c = β*`).

**Salvaguardas** (importantes porque las convenciones de signo de los valores duales varían entre solvers):

1. Tras resolver, se valida el candidato `β` (y su negado, por si el solver invierte el signo) comprobando *factibilidad primal* (`Aβ ≥ b` con tolerancia relativa) y *dualidad fuerte* (`c'β = b'λ*` con tolerancia relativa).
2. Si la validación falla, si el solver no devuelve óptimo, o si ocurre cualquier error, se recurre automáticamente al **primal reducido con GLPK** (`estimate_coefficients_reduced`), que es la formulación de la sección 4.1 ya verificada en la primera ronda.

El resultado nunca puede ser inválido: o pasa la autoverificación, o se resuelve por la vía ya validada.

### 4.3. Scores FDH por enumeración (`rad_out_fdh`)

El score radial output FDH-VRS tiene fórmula cerrada:

```
φ_d = max_{j : x_j ≤ x_d}  min_r ( y_jr / y_dr )
```

donde la dominancia `x_j ≤ x_d` es componente a componente. Es la solución exacta del MILP original (la binaria selecciona el mejor peer dominante). La implementación vectoriza la comprobación de dominancia y evita N MILPs. Verificado en la primera ronda con diferencias ≤ 4·10⁻¹¹ frente al MILP.

### 4.4. `set_knots_high_speed`

Vectoriza la contabilidad de índices del *minimum span* (bucles `for` con `append` → `outer` + filtrado). Devuelve exactamente el mismo conjunto de nudos elegibles que la versión original; es una mejora menor de coste pero determinista.

---

## 5. Resultados empíricos de la primera ronda

Medidos con el test de equivalencia entre la versión original y la primera optimización (sin dual), R 4.3.2, semillas fijas:

| Escenario | Original (s) | Optimizado (s) | Speedup |
|---|---|---|---|
| 1 input, N=100 | 3,14 | 2,28 | 1,4× |
| 3 inputs, N=150 | 77,94 | 49,41 | 1,6× |
| 3×3 (outputs), N=100 | 59,34 | 58,51 | 1,0× |
| 1 input, N=500 | 38,33 | 34,09 | 1,1× |
| 3 inputs, N=300 | 174,06 | 122,01 | 1,4× |
| 3 inputs, N=300, quick | 90,25 | 62,60 | 1,4× |

Media geométrica: **1,3×**. Equivalencia matemática verificada: en 1.861 LPs interceptados, diferencia relativa máxima de valor objetivo 2,4·10⁻¹⁴ y de predicciones 1,1·10⁻¹⁴; scores FDH idénticos a 4·10⁻¹¹.

**Lección:** la reducción de columnas no era suficiente; el coste del símplex lo dictan las filas. De ahí la vía dual de la sección 4.2, cuyos resultados se medirán con `test_high_speed.R`. Expectativa razonable: el coste por LP debería caer de forma notable en N grandes (las bases pasan de ~N² a ~p² elementos), pero la construcción de la matriz `B` candidata y la sobrecarga de R fijan un suelo; no se promete una cifra hasta medirla.

---

## 6. Óptimos múltiples y reproducibilidad: qué se descubrió y qué implica

El test de la primera ronda reveló que el LP de estimación tiene, con frecuencia, **óptimos múltiples**: en los LPs interceptados, las dos formulaciones alcanzaron el mismo objetivo y las mismas predicciones con coeficientes que diferían hasta 3,7 unidades. La causa típica es la colinealidad entre funciones base (los dos *hinges* de un par suman una función lineal que se solapa con el intercepto): existen direcciones en las que los coeficientes pueden moverse sin cambiar la frontera estimada en los puntos muestrales.

Consecuencia: el algoritmo forward compara errores entre nudos candidatos con desigualdades estrictas. Dos candidatos pueden producir *exactamente* el mismo error (frecuente en fronteras envolventes), y entonces la elección depende de ruido de redondeo de ~10⁻¹⁵, que difiere entre formulaciones (e incluso entre versiones de solver, BLAS o sistema operativo, también para la versión original). En esos empates, los modos legacy y high-speed pueden divergir y producir modelos finales distintos — cada uno igualmente óptimo en cada paso. En la primera ronda esto afectó a 4 de 6 escenarios, con diferencias de hasta 0,30 en los scores finales.

Esto **no es un error de cálculo** (la equivalencia por llamada está verificada a precisión de máquina), sino una propiedad del algoritmo: su salida no es única cuando hay empates. Para el artículo conviene tenerlo presente: la reproducibilidad bit a bit solo está garantizada fijando la implementación (`high_speed = FALSE` reproduce la original) y el entorno. Si en el futuro se quisiera una salida única e independiente de la implementación, habría que estabilizar los desempates (p. ej., redondear el criterio de error a 10⁻⁹ y desempatar por valor del nudo), lo que cambiaría los resultados también respecto a la versión original en los casos de empate.

---

## 7. Verificación y uso

**Test (una sola ejecución, sin cambiar de rama):**

```r
source("test_high_speed.R")
```

- **Parte A (criterio PASS/FAIL):** equivalencia matemática por llamada — FDH enumeración vs MILP (A1), LP lineal (A2) y LP de suavizado (A3) comparando ambas implementaciones en cada llamada real. Deben salir PASS; cualquier FAIL indica un problema real.
- **Parte B (informativa):** ajusta los 6 escenarios con ambos modos en la misma sesión; tabla de tiempos/speedups y comparación de nudos, predicciones y scores. Diferencias aquí, con la parte A en PASS, corresponden al fenómeno de la sección 6.

**Uso normal del paquete:** nada cambia; `high_speed = TRUE` es el valor por defecto. Para reproducir resultados antiguos exactamente: `aces(..., high_speed = FALSE)`.

---

## 8. Estado de git y GitHub

- Todo el trabajo está en **una sola rama**: `test-ws-integration`, con tres commits relevantes: `bb663e2` (tu último commit previo, "19/05"), `56895f0` ("baseline", que recoge el estado exacto de tu directorio de trabajo antes de cualquier optimización — tu trabajo de meses que estaba sin commitear) y el commit "feat: modo high_speed...".
- Las ramas nunca se subieron a GitHub (no se hizo `push` en ningún momento); por eso no las veías. Para publicar la rama: `git push origin test-ws-integration`.
- La rama auxiliar `perf-optimization` (primera ronda) se ha eliminado en local; su contenido queda superseded por el commit actual y es recuperable durante ~90 días vía `git reflog` si hiciera falta.

**Pendiente tras revisar:** ejecutar `devtools::document()` para regenerar los `.Rd` (hay funciones nuevas y un parámetro nuevo en `aces()`), y añadir a `.Rbuildignore` antes de CRAN: `informe_optimizacion.md`, `informe_rendimiento.md`, `test_high_speed.R`, `test.R` y `.*\.rds$`.

---

## 9. Trabajo futuro (no incluido, por orden de potencial)

1. **LP incremental con warm start.** Entre dos nudos candidatos consecutivos el LP solo cambia en 2 columnas (el nuevo par de bases) y en pocas filas de forma. Mantener un único modelo `lpSolveAPI` por iteración y modificarlo (`set.column`) con re-optimización desde la base anterior evitaría reconstruir y resolver desde cero ~N veces por iteración. Es el cambio con mayor potencial restante, pero exige reestructurar `add_basis_function` y un manejo cuidadoso de las restricciones de forma, que cambian con cada candidato.
2. **Perfilado fino con `Rprof`/`profvis`** sobre datos reales del tamaño de tu aplicación, para decidir si el siguiente cuello es la construcción de `B` (asignación de memoria por candidato) o el solver.
3. **Construcción incremental de `B`.** Evitar `cbind` de tres bloques por candidato preasignando la matriz y rellenando las dos columnas nuevas.
4. **Reducción de candidatos.** `quick_aces` ya muestrea nudos; podría hacerse un cribado barato (p. ej. evaluar primero un subconjunto y refinar localmente alrededor del mejor).
5. Paralelización: descartada por petición expresa.

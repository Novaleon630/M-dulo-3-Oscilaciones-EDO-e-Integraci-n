# Módulo 3: Oscilaciones por viento (+ apoyo al Módulo 1)

Parte del Integrante 2 del proyecto **SIMTELE-RISK** (Métodos Numéricos). Simulo cómo vibra una estructura bajo la acción del viento, resuelvo la ecuación con Runge-Kutta, integro la magnitud acumulada `x(t)²` y compruebo todo contra soluciones exactas. Además doy apoyo al Módulo 1 (detención de cabinas) con RK2, RK4 e integración.

Todo está en un solo archivo, `modulo3_completo.py`. Los métodos numéricos están programados desde cero, sin NumPy ni SciPy. Solo se usan `flask` (API) y `matplotlib` (gráficos).

---

## 1. Cómo ejecutarlo

```bash
pip install flask matplotlib
python modulo3_completo.py
```

| Comando | Qué hace |
|---|---|
| `python modulo3_completo.py` | Consola, gráficos y API |
| `python modulo3_completo.py --sin-servidor` | Consola y gráficos (sin Flask) |
| `python modulo3_completo.py --sin-graficos` | Consola y API |
| `python modulo3_completo.py --mostrar` | Además abre las ventanas de los gráficos |
| `python modulo3_completo.py --json` | Guarda `resultados_modulo3.json` |

Los gráficos se guardan en la carpeta `graficos_modulo3/`. La API queda en `http://localhost:5001` y se detiene con `Ctrl+C`. Sin `flask` o sin `matplotlib` el archivo sigue funcionando; solo avisa qué falta.

---

## 2. El modelo físico

Una estructura de masa `m`, amortiguamiento `c` y rigidez `k` excitada por una fuerza de viento:

```
m·x''(t) + c·x'(t) + k·x(t) = F(t),     F(t) = F0·sin(wf·t)
```

Para usar Runge-Kutta paso la ecuación de segundo orden a un sistema de primer orden con `y = [x, v]`:

```
x' = v
v' = (F(t) − c·v − k·x) / m
```

Eso es lo que hace `f_oscilador(t, y, m, c, k, F)`. Con `m = 800` y `k = 4000`, la frecuencia natural es `w0 = √(k/m) ≈ 2.236 rad/s`.

---

## 3. Métodos numéricos

**Ecuaciones diferenciales** (firma común `(f, t0, y0, tf, h, params)`):

- `rk2_heun`: Runge-Kutta de 2.º orden, variante de Heun (promedio de la pendiente inicial y la final).
- `rk4`: Runge-Kutta clásico de 4.º orden.

**Integración**, en dos versiones: sobre una función `f(x)` y sobre datos equiespaciados (la solución numérica):

| Sobre una función | Sobre datos | Nota |
|---|---|---|
| `trapecio_compuesto` | `trapecio_datos` | |
| `simpson13_compuesto` | `simpson_datos` | `n` debe ser par |
| `romberg`, `romberg_tabla` | `romberg_datos` | Sobre datos exige `n = 2^m` |

**Errores** (`error_abs`, `error_rel`): el error relativo devuelve `None` si la referencia es cero o casi cero, para no dividir entre cero.

---

## 4. Soluciones exactas (referencia)

- `solucion_analitica`: vibración libre subamortiguada, `x(t) = e^(−ζ·w0·t)·[C1·cos(wd·t) + C2·sin(wd·t)]`.
- `solucion_analitica_forzada`: solución exacta con viento `F0·sin(wf·t)`, igual a la solución particular más la homogénea. Devuelve `(x, v)`. Con `F0 = 0` coincide con la libre.

Requieren `ζ < 1`. Con mis parámetros `ζ ≈ 0.0335` (normal y estrés) y `ζ ≈ 0.00014` (mal condicionado).

Gracias a la solución forzada, **los tres escenarios se comparan contra una solución exacta**, sin usar una referencia numérica.

---

## 5. Los tres escenarios

| | NORMAL | ESTRÉS | MAL CONDICIONADO |
|---|---|---|---|
| `c` | 120 | 120 | 0.5 |
| Viento `F(t)` | 0 | `500·sin(2t)` | `500·sin(0.99·w0·t)` |
| Duración `T` | 10.24 s | 16 s | 40.96 s |
| Subintervalos `n` | 1024 (h = 0.01) | 64 (h = 0.25) | 2048 (h = 0.02) |

Cada `n` es potencia de 2 para poder aplicar Romberg a la solución numérica.

**Normal.** Vibración libre con paso fino. Es el caso de control.

**Estrés.** Viento cerca de la resonancia (2 rad/s frente a 2.236), **paso 25 veces mayor** y más tiempo. Es una exigencia mayor del modelo: RK2 pierde precisión (error y desfase visibles en el gráfico) y RK4 aguanta.

**Mal condicionado.** Es un **subexperimento de sensibilidad**, no un caso de paso grande o falta de convergencia. Con amortiguamiento casi nulo y el viento en `0.99·w0`, la amplitud permanente `X(w) = F0 / √((k − m·w²)² + (c·w)²)` es extremadamente sensible a `wf`: un **1 % de cambio en `wf`** hace pasar la amplitud de **6.28 m a 363.7 m** (se multiplica por 58). La función `sensibilidad_amplitud` devuelve ese factor y el gráfico `m3_mal_condicionado_sensibilidad.png` lo muestra.

Esta distinción importa: el mal condicionamiento es una propiedad del **problema** (cuánto cambia la salida ante un cambio pequeño en los datos), no una falla del método, y no se confunde con inestabilidad, paso grande o riesgo físico.

### Resultados

| Escenario | Error máx. RK2 | Error máx. RK4 | Integral de x² (referencia) | Error de Simpson sobre RK4 |
|---|---|---|---|---|
| Normal | 8.6e-4 | 2.1e-8 | 2.653990 | 1.5e-8 |
| Estrés | 4.8e-1 | 8.3e-3 | 3.353760 | 5.8e-3 |
| Mal condicionado | 8.1e-2 | 8.1e-6 | 212.810845 | 5.3e-6 |

El error máximo se mide sobre `x(t)` en todos los instantes de la malla, contra la solución exacta.

---

## 6. Integración de la magnitud acumulada

La magnitud es `∫₀ᵀ x(t)² dt` (una medida de la energía de la vibración).

- **Valor de referencia:** Simpson 1/3 con `n = 2¹⁷` sobre la solución **exacta**. Con ese `n` el error del método es mucho menor que cualquier diferencia que se reporta.
- **Sobre la simulación:** trapecio, Simpson y Romberg aplicados a `x²` calculado por RK4 (`integral`) y por RK2 (`integral_rk2`), con su error absoluto y relativo (`integral_errores`).
- **Tabla de refinamiento** (`refinamiento`): para `n = 2, 4, …, 4096`, el valor y el error de cada método sobre la función exacta. En el escenario normal, con `n = 4096`, los errores son 2.6e-7 (trapecio), 1.7e-12 (Simpson) y 3.8e-14 (Romberg). Se ve el orden de cada método.

En el escenario de estrés el error de la integral sobre RK4 (≈ 6e-3) lo domina el error de la EDO, no el de la integración; por eso Romberg no mejora a Simpson ahí.

---

## 7. Apoyo al Módulo 1 (detención de cabinas)

`modulo1_apoyo` resuelve con RK2 y RK4 el frenado `x' = v`, `v' = −a0·(1 + β·t)`. Con `β = 0` la desaceleración es constante; con `β > 0` es variable. Luego calcula la posición **integrando `v(t)`** con trapecio, Simpson y Romberg y la compara con la solución exacta.

| Caso | `t_det` | Distancia de frenado | Distancia total (con reacción 0.5 s) |
|---|---|---|---|
| Constante (`v0 = 12`, `a0 = 2.5`) | 4.800 s | 28.800 m | 34.800 m |
| Variable (`β = 0.2`) | 3.544 s | 23.119 m | 29.119 m |

En el caso constante coincide con `t = v0/a` y `d = v0²/(2a)`. El tiempo de detención se toma analítico; el Integrante 1 puede reemplazarlo por la raíz que calculen sus métodos (bisección, Newton o secante).

---

## 8. Pruebas con solución analítica

`pruebas_analiticas()` ejecuta y reporta `OK` o `FALLA`:

1. `∫₀² x² dx = 8/3` con trapecio (n = 8, error 0.0208 como predice la teoría), Simpson y Romberg (exactos).
2. Orden de convergencia con `y' = −y`: al reducir `h` a la mitad, el error de RK2 baja ≈ 4 veces y el de RK4 ≈ 16 veces.
3. RK4 contra la solución libre exacta (error < 1e-6).
4. Solución analítica forzada contra RK4 con `h = 0.001` (coinciden a 3e-13).
5. La segunda columna de Romberg coincide con Simpson.

También están disponibles en `GET /api/pruebas`.

---

## 9. API y formato de salida

| Ruta | Respuesta |
|---|---|
| `GET /api/modulo3/<tipo>` | `NORMAL`, `ESTRES` o `MAL_CONDICIONADO` (no distingue mayúsculas). Otro valor: `400` |
| `GET /api/modulo1_apoyo/<tipo>` | `constante` o `variable` |
| `GET /api/pruebas` | Resultado de las pruebas analíticas |

Los resultados se calculan una sola vez y se guardan en memoria.

Campos principales de `/api/modulo3/<tipo>`:

| Campo | Contenido |
|---|---|
| `escenario`, `parametros` | Nombre, parámetros, `h`, `w0`, `zeta` |
| `rk2`, `rk4` | `ts`, `xs`, `vs`, `errores`, `error_max`, `amplitud_max`, `amplificacion`, `paso`, `convergio`, `tiempo_ms` |
| `analitica` | Solución exacta `x(t)` para graficar |
| `iteraciones` | Filas para la tabla `iteracion_resultado` (hasta 200 por método): `numero`, `aproximacion` `{t, x, v}`, `residual`, `error_abs`, `error_rel` |
| `integral`, `integral_rk2`, `integral_errores`, `integral_referencia` | Integrales de `x²` y sus errores |
| `refinamiento` | Tabla de refinamiento de los tres métodos |
| `sensibilidad` | Solo en `MAL_CONDICIONADO` |

Notas para quien conecte la base de datos:

- `convergio` significa que la simulación fue **estable** (valores finitos y sin crecer sin control), no que haya una tolerancia de convergencia.
- En `iteraciones`, `residual` es el error de la segunda variable de estado (la velocidad), porque en un Runge-Kutta no existe un residual de la ecuación.
- `error_rel` puede ser `null` cuando la referencia es cero.

---

## 10. Gráficos

Se generan en `graficos_modulo3/`, al menos dos por escenario:

- `m3_<escenario>_solucion_error.png`: `x(t)` exacta, RK2 y RK4 superpuestas, y el error absoluto.
- `m3_<escenario>_integracion.png`: error de la integral frente a `n` para trapecio, Simpson y Romberg.
- `m3_mal_condicionado_sensibilidad.png`: amplitud permanente `X(w)` y el efecto de un 1 % en `wf`.
- `m1_apoyo_constante.png` y `m1_apoyo_variable.png`: `v(t)` y posición durante el frenado.

---

## 11. Limitaciones

- El modelo es de un solo grado de libertad y lineal. No representa el comportamiento real de una cabina ni de una estructura; es un modelo educativo con datos sintéticos.
- La sensibilidad del escenario mal condicionado se mide sobre la amplitud en régimen permanente (fórmula exacta). Con `c = 0.5` el sistema tardaría mucho en llegar a ese régimen, así que la simulación de 40.96 s muestra la respuesta transitoria; el factor de amplificación corresponde al problema matemático, no a lo que se ve en esa ventana de tiempo.
- RK2 no se vuelve inestable en el escenario de estrés (el amortiguamiento lo contiene), solo pierde precisión. Con amortiguamiento casi nulo y un paso grande sí puede amplificar la amplitud; eso sería un problema de estabilidad del método, distinto del mal condicionamiento.

---

## 12. Estructura del archivo

```
0) UTILIDADES DE ERROR     error_abs, error_rel
1) MÉTODOS NUMÉRICOS       f_oscilador, rk2_heun, rk4, trapecio, simpson, romberg (+ versiones sobre datos)
2) MODELO Y SOL. EXACTAS   solucion_analitica, solucion_analitica_forzada, sensibilidad_amplitud
3) ESCENARIOS MÓDULO 3     PARAMETROS, ejecutar_escenario, escenario_normal / estres / mal_condicionado
4) APOYO MÓDULO 1          modulo1_apoyo, casos_modulo1
5) PRUEBAS ANALÍTICAS      pruebas_analiticas
6) CONSOLA                 probar_en_consola
7) GRÁFICOS                generar_graficos
8) JSON Y CACHÉ            obtener_todos, guardar_json
9) API                     Flask: /api/modulo3, /api/modulo1_apoyo, /api/pruebas
```

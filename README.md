# Módulo 3: Oscilaciones por viento

En este módulo simulo cómo vibra una estructura (por ejemplo, un puente o una torre) cuando el viento la excita. Todo está en un solo archivo de Python: implementé yo mismo los métodos numéricos desde cero (sin NumPy ni SciPy), armé tres escenarios de prueba y expuse los resultados mediante una pequeña API con Flask para que mi compañero los consuma desde la web.

---

## 1. El modelo físico

Modelo la estructura como un **oscilador amortiguado con fuerza externa**:

```
m·x''(t) + c·x'(t) + k·x(t) = F(t)
```

| Símbolo | Significado | Valor que uso |
|---|---|---|
| `m` | masa | 800 |
| `c` | amortiguamiento | 120 (normal y estrés) / 0.5 (mal condicionado) |
| `k` | rigidez del resorte | 4000 |
| `F(t)` | fuerza del viento | 0 o `500·sin(2t)` según el escenario |
| `x` | desplazamiento | incógnita |
| `v = x'` | velocidad | incógnita |

Los métodos de Runge-Kutta trabajan con ecuaciones de **primer orden**, así que convierto la ecuación de segundo orden en un sistema de dos ecuaciones con el vector `y = [x, v]`:

```
x' = v
v' = (F(t) − c·v − k·x) / m
```

Esto es exactamente lo que hace la función `f_oscilador(t, y, m, c, k, F)`: recibe el estado actual y devuelve sus derivadas `[v, (F(t) - c*v - k*x)/m]`.

---

## 2. Los métodos numéricos (sección 1 del código)

### Para resolver la ecuación diferencial

- **`rk2_heun`**: Runge-Kutta de 2.º orden (método de Heun). En cada paso calculo la pendiente al inicio (`k1`), doy un paso de prueba tipo Euler (`y_pred`), calculo la pendiente al final (`k2`) y avanzo con el **promedio** de las dos pendientes.
- **`rk4`**: Runge-Kutta clásico de 4.º orden. Calculo cuatro pendientes (`k1` a `k4`) y avanzo con el promedio ponderado `(k1 + 2·k2 + 2·k3 + k4)/6`. Es más preciso que Heun con el mismo paso `h`, aunque cuesta más evaluaciones de `f` por paso.

Las dos funciones tienen la misma firma: `(f, t0, y0, tf, h, params)`. Devuelven la lista de tiempos `ts` y la lista de estados `ys`. El número de pasos lo calculo con `n = round((tf - t0)/h)`.

### Para calcular integrales

Más adelante necesito integrar `x(t)²` (una medida de la energía de la vibración), así que implementé tres métodos:

- **`trapecio_compuesto(f, a, b, n)`**: aproxima el área con `n` trapecios.
- **`simpson13_compuesto(f, a, b, n)`**: usa parábolas en lugar de rectas. Exige que `n` sea **par**; si no, lanza `ValueError`.
- **`romberg(f, a, b, nivel_max, tol)`**: parte del trapecio y aplica **extrapolación de Richardson** para ir mejorando la precisión. Se detiene antes si la diferencia entre dos niveles consecutivos de la diagonal es menor que `tol`.

---

## 3. La solución analítica (mi "verdad" de referencia)

`solucion_analitica(t, x0, v0, m, c, k)` calcula la solución exacta del caso **sin fuerza externa y subamortiguado**:

```
ω0 = √(k/m)                 (frecuencia natural)
ζ  = c / (2·√(k·m))         (razón de amortiguamiento)
ωd = ω0·√(1 − ζ²)           (frecuencia amortiguada)

x(t) = e^(−ζ·ω0·t) · [ C1·cos(ωd·t) + C2·sin(ωd·t) ]
C1 = x0,    C2 = (v0 + ζ·ω0·x0) / ωd
```

Si `ζ ≥ 1` el sistema ya no oscila (sobreamortiguado o crítico) y esta fórmula no aplica, por eso lanzo un error en ese caso. Con mis parámetros, `ζ ≈ 0.034` en los escenarios normal y de estrés, y `ζ ≈ 0.00014` en el mal condicionado, así que siempre estoy en el caso subamortiguado.

---

## 4. Los tres escenarios (sección 2 del código)

Cada escenario es una función que corre RK2 y RK4, mide el error, calcula las integrales de `x²` y devuelve todo en un diccionario.

### `escenario_normal()`: el caso "sano"

- Condiciones: `x0 = 1`, `v0 = 0`, sin viento (`F = 0`), `h = 0.01`, de `t = 0` a `10`.
- Como aquí existe solución exacta, el **error máximo** de RK2 y RK4 lo calculo contra `solucion_analitica` en cada punto.
- Integro `x²` con Trapecio y Simpson (`n = 1000`) y con Romberg (nivel 6).
- Es el escenario donde espero ver que **RK4 es muchísimo más preciso que RK2**.

### `escenario_estres()`: viento fuerte cerca de la resonancia

- Condiciones: parto del reposo (`x0 = 0`, `v0 = 0`) y aplico `F(t) = 500·sin(2t)`. La frecuencia del viento (2 rad/s) está muy cerca de la frecuencia natural de la estructura (≈ 2.236 rad/s), lo que provoca una amplificación fuerte.
- Uso un paso **grande** (`h = 0.25`) a propósito, para ver cómo se degradan los métodos.
- Con fuerza externa no tengo solución analítica (por eso `"analitica": []`). Como referencia calculo una solución **RK4 con paso muy fino** (`h = 0.001`) y comparo contra ella en los mismos tiempos.
- Integración con `n = 200` y Romberg nivel 5.

### `escenario_mal_condicionado()`: casi sin amortiguamiento y paso grande

- Condiciones: `c = 0.5` (amortiguamiento casi nulo), `h = 0.30` (paso grande), sin viento.
- Es un caso difícil: la oscilación casi no se disipa y el paso es grande respecto al período, por lo que los errores numéricos se acumulan y se notan mucho.
- El error se mide contra la referencia RK4 con `h = 0.001`; además devuelvo la solución analítica para graficar.
- Integración con `n = 300` y Romberg nivel 5.

**Sobre la referencia fina:** para comparar, tomo el valor de la solución de referencia en el índice `round(t/0.001)` de cada tiempo `t` de la malla gruesa.

---

## 5. Prueba por consola (sección 3 del código)

`probar_en_consola()` ejecuta los tres escenarios e imprime, para cada uno:

- Error máximo de RK2
- Error máximo de RK4
- Integral de `x²` por Trapecio, Simpson y Romberg

La uso para verificar rápido que los números tienen sentido antes de conectar la web.

---

## 6. La API (sección 4 del código)

Creo una app Flask con una sola ruta:

```
GET /api/modulo3/<tipo>
```

| `<tipo>` | Escenario |
|---|---|
| `NORMAL` | `escenario_normal()` |
| `ESTRES` | `escenario_estres()` |
| `MAL_CONDICIONADO` | `escenario_mal_condicionado()` |

El tipo no distingue mayúsculas (hago `tipo.upper()`). Si mandan otro valor, respondo `400` con `{"error": "Escenario no válido"}`.

### Formato de la respuesta (JSON)

```json
{
  "escenario": "NORMAL",
  "rk2": { "ts": [...], "xs": [...], "error_max": 0.0 },
  "rk4": { "ts": [...], "xs": [...], "error_max": 0.0 },
  "integral": { "trapecio": 0.0, "simpson": 0.0, "romberg": 0.0 },
  "analitica": [...]
}
```

- `ts`: tiempos; `xs`: desplazamiento calculado en cada tiempo (lo que se grafica).
- `error_max`: error máximo absoluto respecto a la referencia del escenario.
- `integral`: valor aproximado de la integral de `x²` en `[0, 10]` con cada método.
- `analitica`: solución exacta para graficar (vacía en el escenario de estrés).

---

## 7. Cómo ejecutarlo

Necesito Python 3 y Flask:

```bash
pip install flask
python modulo3.py
```

(Sustituye `modulo3.py` por el nombre real del archivo.)

Al ejecutarlo pasa esto, en orden:

1. Se imprimen en consola los resultados de los tres escenarios.
2. Se levanta el servidor en `http://localhost:5001`.

Para probar la API desde el navegador:

- http://localhost:5001/api/modulo3/NORMAL
- http://localhost:5001/api/modulo3/ESTRES
- http://localhost:5001/api/modulo3/MAL_CONDICIONADO

---

## 8. Cosas que tengo presentes

- **Reinicio en modo debug:** como uso `debug=True`, Flask reinicia el proceso y la salida de consola puede aparecer **dos veces**. Es normal.
- **CORS:** si la página web se sirve desde otro puerto u origen, el navegador puede bloquear las peticiones. En ese caso instalaría `flask-cors` y haría `CORS(app)`.
- **Integral en el escenario de estrés:** ahí integro `x²` usando la solución analítica **libre** (la de `x0 = 1`, sin viento), no la solución forzada. Es un punto a revisar si quiero que la integral represente la energía real bajo el viento; la alternativa sería integrar sobre la referencia RK4 fina.
- **Costo de cómputo:** cada petición recalcula todo (incluida la referencia con `h = 0.001`). Si se vuelve lento, puedo guardar los resultados en caché al arrancar.

---

## 9. Estructura del archivo

```
1) MÉTODOS NUMÉRICOS     f_oscilador, rk2_heun, rk4, trapecio, simpson, romberg
2) MODELO Y ESCENARIOS   solucion_analitica, escenario_normal / estres / mal_condicionado
3) CONSOLA               probar_en_consola
4) API                   Flask, ruta /api/modulo3/<tipo>
```

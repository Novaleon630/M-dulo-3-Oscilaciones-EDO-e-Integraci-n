import math
from flask import Flask, jsonify

# ============================================================
# 1) MÉTODOS NUMÉRICOS (Desde cero)
# ============================================================

def f_oscilador(t, y, m, c, k, F):
    """Ecuación diferencial del oscilador: y = [x, v]"""
    x, v = y
    return [v, (F(t) - c * v - k * x) / m]

def rk2_heun(f, t0, y0, tf, h, params):
    """Runge-Kutta de 2.º orden (Heun)"""
    ts, ys = [t0], [list(y0)]
    t, y = t0, list(y0)
    n = int(round((tf - t0) / h))
    for _ in range(n):
        k1 = f(t, y, *params)
        y_pred = [y[i] + h * k1[i] for i in range(len(y))]
        k2 = f(t + h, y_pred, *params)
        y = [y[i] + (h / 2.0) * (k1[i] + k2[i]) for i in range(len(y))]
        t += h
        ts.append(t)
        ys.append(list(y))
    return ts, ys

def rk4(f, t0, y0, tf, h, params):
    """Runge-Kutta de 4.º orden clásico"""
    ts, ys = [t0], [list(y0)]
    t, y = t0, list(y0)
    n = int(round((tf - t0) / h))
    for _ in range(n):
        k1 = f(t, y, *params)
        y2 = [y[i] + 0.5 * h * k1[i] for i in range(len(y))]
        k2 = f(t + 0.5 * h, y2, *params)
        y3 = [y[i] + 0.5 * h * k2[i] for i in range(len(y))]
        k3 = f(t + 0.5 * h, y3, *params)
        y4 = [y[i] + h * k3[i] for i in range(len(y))]
        k4 = f(t + h, y4, *params)
        y = [y[i] + (h / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(len(y))]
        t += h
        ts.append(t)
        ys.append(list(y))
    return ts, ys

def trapecio_compuesto(f, a, b, n):
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        s += 2 * f(a + i * h)
    return (h / 2.0) * s

def simpson13_compuesto(f, a, b, n):
    if n % 2 != 0:
        raise ValueError("n debe ser par")
    h = (b - a) / n
    s = f(a) + f(b)
    for i in range(1, n):
        s += (4 if i % 2 == 1 else 2) * f(a + i * h)
    return (h / 3.0) * s

def romberg(f, a, b, nivel_max=5, tol=1e-12):
    R = [[0.0] * (nivel_max + 1) for _ in range(nivel_max + 1)]
    h = b - a
    R[0][0] = 0.5 * h * (f(a) + f(b))
    for i in range(1, nivel_max + 1):
        h /= 2.0
        suma = sum(f(a + k * h) for k in range(1, 2**i, 2))
        R[i][0] = 0.5 * R[i-1][0] + h * suma
        for j in range(1, i + 1):
            R[i][j] = R[i][j-1] + (R[i][j-1] - R[i-1][j-1]) / (4**j - 1)
        if abs(R[i][i] - R[i-1][i-1]) < tol:
            return R[i][i]
    return R[nivel_max][nivel_max]


# ============================================================
# 2) MODELO Y ESCENARIOS
# ============================================================

def solucion_analitica(t, x0, v0, m, c, k):
    w0 = math.sqrt(k / m)
    zeta = c / (2.0 * math.sqrt(k * m))
    if zeta >= 1.0:
        raise ValueError("Requiere zeta < 1")
    wd = w0 * math.sqrt(1.0 - zeta**2)
    C1 = x0
    C2 = (v0 + zeta * w0 * x0) / wd
    return math.exp(-zeta * w0 * t) * (C1 * math.cos(wd * t) + C2 * math.sin(wd * t))

def escenario_normal():
    m, c, k = 800.0, 120.0, 4000.0
    x0, v0, t0, tf, h = 1.0, 0.0, 0.0, 10.0, 0.01
    F = lambda t: 0.0
    params = (m, c, k, F)
    
    ts2, ys2 = rk2_heun(f_oscilador, t0, [x0, v0], tf, h, params)
    ts4, ys4 = rk4(f_oscilador, t0, [x0, v0], tf, h, params)
    
    err2 = max(abs(ys2[i][0] - solucion_analitica(ts2[i], x0, v0, m, c, k)) for i in range(len(ts2)))
    err4 = max(abs(ys4[i][0] - solucion_analitica(ts4[i], x0, v0, m, c, k)) for i in range(len(ts4)))
    
    f_int = lambda t: solucion_analitica(t, x0, v0, m, c, k)**2
    return {
        "escenario": "NORMAL",
        "rk2": {"ts": ts2, "xs": [y[0] for y in ys2], "error_max": err2},
        "rk4": {"ts": ts4, "xs": [y[0] for y in ys4], "error_max": err4},
        "integral": {
            "trapecio": trapecio_compuesto(f_int, t0, tf, 1000),
            "simpson": simpson13_compuesto(f_int, t0, tf, 1000),
            "romberg": romberg(f_int, t0, tf, 6)
        },
        "analitica": [solucion_analitica(t, x0, v0, m, c, k) for t in ts4]
    }

def escenario_estres():
    m, c, k = 800.0, 120.0, 4000.0
    x0, v0, t0, tf, h = 0.0, 0.0, 0.0, 10.0, 0.25
    F = lambda t: 500.0 * math.sin(2.0 * t)
    params = (m, c, k, F)
    
    ts2, ys2 = rk2_heun(f_oscilador, t0, [x0, v0], tf, h, params)
    ts4, ys4 = rk4(f_oscilador, t0, [x0, v0], tf, h, params)
    ts_ref, ys_ref = rk4(f_oscilador, t0, [x0, v0], tf, 0.001, params)
    
    ref_vals = [ys_ref[min(int(round(t/0.001)), len(ys_ref)-1)][0] for t in ts4]
    err2 = max(abs(ys2[i][0] - ref_vals[i]) for i in range(len(ts2)))
    err4 = max(abs(ys4[i][0] - ref_vals[i]) for i in range(len(ts4)))
    
    f_int = lambda t: solucion_analitica(t, 1.0, 0.0, m, c, k)**2
    return {
        "escenario": "ESTRES",
        "rk2": {"ts": ts2, "xs": [y[0] for y in ys2], "error_max": err2},
        "rk4": {"ts": ts4, "xs": [y[0] for y in ys4], "error_max": err4},
        "integral": {
            "trapecio": trapecio_compuesto(f_int, t0, tf, 200),
            "simpson": simpson13_compuesto(f_int, t0, tf, 200),
            "romberg": romberg(f_int, t0, tf, 5)
        },
        "analitica": []
    }

def escenario_mal_condicionado():
    m, c, k = 800.0, 0.5, 4000.0
    x0, v0, t0, tf, h = 1.0, 0.0, 0.0, 10.0, 0.30
    F = lambda t: 0.0
    params = (m, c, k, F)
    
    ts2, ys2 = rk2_heun(f_oscilador, t0, [x0, v0], tf, h, params)
    ts4, ys4 = rk4(f_oscilador, t0, [x0, v0], tf, h, params)
    ts_ref, ys_ref = rk4(f_oscilador, t0, [x0, v0], tf, 0.001, params)
    
    ref_vals = [ys_ref[min(int(round(t/0.001)), len(ys_ref)-1)][0] for t in ts4]
    err2 = max(abs(ys2[i][0] - ref_vals[i]) for i in range(len(ts2)))
    err4 = max(abs(ys4[i][0] - ref_vals[i]) for i in range(len(ts4)))
    
    f_int = lambda t: solucion_analitica(t, x0, v0, m, c, k)**2
    return {
        "escenario": "MAL_CONDICIONADO",
        "rk2": {"ts": ts2, "xs": [y[0] for y in ys2], "error_max": err2},
        "rk4": {"ts": ts4, "xs": [y[0] for y in ys4], "error_max": err4},
        "integral": {
            "trapecio": trapecio_compuesto(f_int, t0, tf, 300),
            "simpson": simpson13_compuesto(f_int, t0, tf, 300),
            "romberg": romberg(f_int, t0, tf, 5)
        },
        "analitica": [solucion_analitica(t, x0, v0, m, c, k) for t in ts4]
    }


# ============================================================
# 3) IMPRESIÓN EN CONSOLA (Para que pruebes)
# ============================================================

def probar_en_consola():
    print("=" * 60)
    print("MÓDULO 3: OSCILACIONES POR VIENTO")
    print("=" * 60)
    
    escenarios = [
        ("NORMAL", escenario_normal()),
        ("ESTRÉS", escenario_estres()),
        ("MAL CONDICIONADO", escenario_mal_condicionado())
    ]
    
    for nombre, res in escenarios:
        print(f"\n--- Escenario {nombre} ---")
        print(f"  Error máximo RK2: {res['rk2']['error_max']:.6e}")
        print(f"  Error máximo RK4: {res['rk4']['error_max']:.6e}")
        print(f"  Integral x^2 (Trapecio): {res['integral']['trapecio']:.6f}")
        print(f"  Integral x^2 (Simpson) : {res['integral']['simpson']:.6f}")
        print(f"  Integral x^2 (Romberg) : {res['integral']['romberg']:.6f}")


# ============================================================
# 4) API PARA TU COMPAÑERO (Flask)
# ============================================================

app = Flask(__name__)

@app.route("/api/modulo3/<tipo>")
def api_modulo3(tipo):
    tipo = tipo.upper()
    if tipo == "NORMAL":
        return jsonify(escenario_normal())
    elif tipo == "ESTRES":
        return jsonify(escenario_estres())
    elif tipo == "MAL_CONDICIONADO":
        return jsonify(escenario_mal_condicionado())
    else:
        return jsonify({"error": "Escenario no válido"}), 400


if __name__ == "__main__":
    # Al ejecutar este archivo, primero muestra los resultados en consola
    probar_en_consola()
    
    # Luego levanta el servidor web para que tu compañero lo pruebe
    print("\n" + "=" * 60)
    print("Servidor API corriendo en http://localhost:5001")
    print("Presiona Ctrl+C para detenerlo")
    print("=" * 60)
    app.run(debug=True, port=5001)
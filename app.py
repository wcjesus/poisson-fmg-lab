#!/usr/bin/env python3
"""
Interface web (Flask) para o LABORATÓRIO MULTIGRID (executável `lab`).

Estrutura esperada:

    pasta/
      app.py
      templates/index.html
      lab/                 <- código-fonte C (Makefile, *.c, *.h ...)

Outra pasta de fontes:  LAB_DIR=/caminho/lab python3 app.py

O `lab` é interativo (lê A, B, C, pesos, tolerâncias... do teclado). O app o
executa por stdin, sem alterar seus fontes nem seu config.mk. Para cada valor de
EXP (macro de compilação) ele compila um executável próprio em _cache/ e roda
cada requisição numa pasta temporária.

Observação importante sobre o programa original: cada execução resolve a malha N
E a malha 2N (para estimar a ordem p); o gráfico da solução é o da malha 2N.
"""
import base64
import glob
import io
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from flask import Flask, jsonify, render_template, request

app = Flask(__name__)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.environ.get("LAB_DIR") or next(
    (d for d in (os.path.join(BASE_DIR, "lab"), BASE_DIR)
     if os.path.isfile(os.path.join(d, "Makefile"))),
    os.path.join(BASE_DIR, "lab"),
)
CACHE_DIR = os.path.join(BASE_DIR, "_cache")

TIMEOUT_EXEC = 120          # segundos por execução
TIMEOUT_BUILD = 240
_sem_exec = threading.Semaphore(4)
_trava_build = {}
_trava_global = threading.Lock()

KMAX = 5                    # N = 16*2^k, k <= 5  -> N <= 512 (gráfico: 2N = 1024)
KMAX_COMPARAR = 4           # N <= 256 no modo "50 combinações"

NOMES_A = {1: "Jacobi", 2: "Gauss-Seidel"}
NOMES_B = {1: "LU", 2: "LDL", 3: "Gauss-Jordan", 4: "Jacobi", 5: "Gauss-Seidel"}
NOMES_C = {1: "Jacobi", 2: "Gauss-Seidel", 3: "LU", 4: "LDL", 5: "Gauss-Jordan"}


# ----------------------------------------------------------------------
# Valores padrão (lidos de configuracao.h quando possível)
# ----------------------------------------------------------------------
def _define(nome, padrao):
    try:
        with open(os.path.join(SRC_DIR, "configuracao.h"), errors="replace") as f:
            for linha in f:
                m = re.match(rf"\s*#define\s+{nome}\s+(.*?)\s*(//.*)?$", linha)
                if m:
                    expr = m.group(1).strip()
                    if re.fullmatch(r"[0-9eE.+\-*/() ]+", expr):
                        return float(eval(expr, {"__builtins__": {}}))
    except (OSError, SyntaxError, ZeroDivisionError):
        pass
    return padrao


def _exp10(x, padrao):
    return int(round(-math.log10(x))) if x and x > 0 else padrao


PADRAO = {
    "K": _define("K", 1.0),
    "EXP": int(_define("EXP", 1)),
    "OMEGA": _define("OMEGA", 0.8),            # Jacobi
    "OMEGAGS": _define("OMEGAGS", 1.0),        # Gauss-Seidel
    "NU1": int(_define("NU1", 1)),
    "NU2": int(_define("NU2", 1)),
    "TOL": _exp10(_define("T", 1e-7), 7),
    "TOLB": _exp10(_define("TOLERANCIA_B", 1e-10), 10),
    "TOLC": _exp10(_define("TOLERANCIA_C", 1e-2), 2),
}
FAIXAS = {"K": (1.0, 8.0), "EXP": (1, 3), "OMEGA": (0.1, 1.95),
          "NU": (1, 6), "TOL": (1, 14)}


# ----------------------------------------------------------------------
# Validação
# ----------------------------------------------------------------------
def _int(nome, valor, lo, hi):
    try:
        n = int(valor)
    except (TypeError, ValueError):
        raise ValueError(f"{nome} deve ser inteiro.")
    if not lo <= n <= hi:
        raise ValueError(f"{nome} deve estar entre {lo} e {hi}.")
    return n


def _flt(nome, valor, lo, hi):
    try:
        x = float(valor)
    except (TypeError, ValueError):
        raise ValueError(f"{nome} deve ser numérico.")
    if not (lo <= x <= hi) or math.isnan(x):
        raise ValueError(f"{nome} deve estar entre {lo:g} e {hi:g}.")
    return x


def _padrao_omega(metodo_iterativo_gs):
    return PADRAO["OMEGAGS"] if metodo_iterativo_gs else PADRAO["OMEGA"]


def _parametros(f, comparar=False):
    kmax = KMAX_COMPARAR if comparar else KMAX
    p = {}
    p["k"] = _int("malha", f.get("k", 3), 0, kmax)
    p["N"] = 16 * 2 ** p["k"]
    p["A"] = _int("A", f.get("A", 2), 1, 2)
    p["B"] = _int("B", f.get("B", 1), 1, 5)
    p["C"] = _int("C", f.get("C", 1), 1, 5)
    om = FAIXAS["OMEGA"]
    p["omA"] = _flt("ω_A", f.get("omA", _padrao_omega(p["A"] == 2)), *om)
    p["omB"] = _flt("ω_B", f.get("omB", _padrao_omega(p["B"] == 5)), *om)
    p["omC"] = _flt("ω_C", f.get("omC", _padrao_omega(p["C"] == 2)), *om)
    p["nu1"] = _int("NU1", f.get("nu1", PADRAO["NU1"]), *FAIXAS["NU"])
    p["nu2"] = _int("NU2", f.get("nu2", PADRAO["NU2"]), *FAIXAS["NU"])
    p["tol"] = 10.0 ** -_int("tolerância global", f.get("tol", PADRAO["TOL"]), *FAIXAS["TOL"])
    p["tolB"] = 10.0 ** -_int("tolerância de B", f.get("tolB", PADRAO["TOLB"]), *FAIXAS["TOL"])
    p["tolC"] = 10.0 ** -_int("tolerância de C", f.get("tolC", PADRAO["TOLC"]), *FAIXAS["TOL"])
    p["K"] = _flt("K", f.get("K", PADRAO["K"]), *FAIXAS["K"])
    p["EXP"] = _int("EXP", f.get("EXP", PADRAO["EXP"]), *FAIXAS["EXP"])
    return p


# ----------------------------------------------------------------------
# Compilação (uma por valor de EXP)
# ----------------------------------------------------------------------
def _garante_executaveis(exp):
    chave = f"EXP{exp}"
    with _trava_global:
        trava = _trava_build.setdefault(chave, threading.Lock())
    pasta = os.path.join(CACHE_DIR, chave)
    with trava:
        os.makedirs(pasta, exist_ok=True)
        with open(os.path.join(SRC_DIR, "configuracao.h"), newline="", errors="replace") as f:
            texto = f.read()
        texto, n = re.subn(r"(?m)^(\s*#define\s+EXP\s+).*?(\r?)$",
                           lambda m: f"{m.group(1)}{exp}{m.group(2)}", texto)
        if n != 1:
            raise RuntimeError("#define EXP não encontrado em configuracao.h.")
        for arq in glob.glob(os.path.join(SRC_DIR, "*.c")) + \
                   glob.glob(os.path.join(SRC_DIR, "*.h")) + \
                   [os.path.join(SRC_DIR, "Makefile")]:
            nome = os.path.basename(arq)
            if nome == "configuracao.h":
                continue
            destino = os.path.join(pasta, nome)
            if not os.path.lexists(destino):
                os.symlink(arq, destino)
        destino = os.path.join(pasta, "configuracao.h")
        if os.path.islink(destino):
            os.remove(destino)
        antigo = open(destino, newline="").read() if os.path.exists(destino) else None
        if antigo != texto:
            with open(destino, "w", newline="") as f:
                f.write(texto)
        with open(os.path.join(pasta, "config.mk"), "w") as f:
            f.write("PRINT = ON\n")
        r = subprocess.run(["make", "-j4", "lab", "labplot"], cwd=pasta,
                           capture_output=True, text=True, timeout=TIMEOUT_BUILD)
        if r.returncode != 0:
            raise RuntimeError("Falha ao compilar:\n" + (r.stderr or r.stdout)[-2000:])
    return os.path.join(pasta, "lab"), os.path.join(pasta, "labplot")


# ----------------------------------------------------------------------
# Execução
# ----------------------------------------------------------------------
def _g(x):
    return repr(float(x))


def monta_entrada(p):
    """Respostas, na ordem em que lab.c as pede (modo 'não usar padrão')."""
    L = ["1", "0", str(p["A"]), _g(p["omA"]), str(p["nu1"]), str(p["nu2"]), str(p["B"])]
    if p["B"] >= 4:
        L += [_g(p["omB"]), _g(p["tolB"])]
    L.append(str(p["C"]))
    if p["C"] <= 2:
        L += [_g(p["omC"]), _g(p["tolC"])]
    L += [_g(p["K"]), _g(p["tol"])]
    return "\n".join(L) + "\n"


def _roda(exe, N, entrada, cwd, timeout=TIMEOUT_EXEC):
    with _sem_exec:
        r = subprocess.run([exe, str(N), "1"], cwd=cwd, input=entrada.encode(),
                           capture_output=True, timeout=timeout)
    return (r.returncode, r.stdout.decode("utf-8", errors="replace"),
            r.stderr.decode("utf-8", errors="replace"))


# ----------------------------------------------------------------------
# Análise da saída
# ----------------------------------------------------------------------
_NUM = r"(-?nan|-?inf|[0-9.eE+\-]+)"


def _f(m):
    return float(m.group(1)) if m else float("nan")


def _i(m):
    return int(m.group(1)) if m else 0


def parse_execucoes(out):
    """Uma entrada por chamada de Multigrid (N e 2N)."""
    blocos = re.split(r"\nMULTIGRID - FMG\n", out)[1:]
    runs = []
    for b in blocos:
        s = lambda pad: re.search(pad, b)
        runs.append({
            "nx": _i(s(r"nx = (\d+)")),
            "ciclos_fmg": _i(s(r"Ciclos V-cycle apos FMG:\s*(\d+)")),
            "iter_C": _i(s(r"Iteracoes em C:\s*(\d+)")),
            "residuo": _f(s(rf"Residuo final:\s*{_NUM}")),
            "erro2": _f(s(rf"Norma 2 do erro:\s*{_NUM}")),
            "erro_inf": _f(s(rf"Norma infinito do erro:\s*{_NUM}")),
            "tempo": _f(s(rf"Tempo de execucao:\s*{_NUM} s")),
            "convergiu": "Convergencia atingida" in b,
            "t_guess": _f(s(rf"Tempo pra gerar o chute:\s*{_NUM} s")),
            "t_vcycle": _f(s(rf"v-cycles \(loop principal\):\s*{_NUM} s")),
            "t_pre": _f(s(rf"Pre-smoothing : total {_NUM} s")),
            "n_pre": _i(s(r"Pre-smoothing : total \S+ s\s+\((\d+) chamadas")),
            "t_pos": _f(s(rf"Pos-smoothing : total {_NUM} s")),
            "n_pos": _i(s(r"Pos-smoothing : total \S+ s\s+\((\d+) chamadas")),
            "t_B": _f(s(rf"B \(coarse\): total {_NUM} s")),
            "n_B_stages": _i(s(r"B \(coarse\): total \S+ s\s+\((\d+) etapas")),
            "n_B_iter": _i(s(r"B iterativo: (\d+) iteracoes")),
            "n_vcycle": _i(s(r"V-cycle: (\d+) ciclos")),
            "hist": [float(v) for _, v in re.findall(
                rf"(?m)^Ciclo (\d+) res[ií]duo: {_NUM}", b)],
        })
    return runs


def parse_ordem(out):
    m = re.search(rf"p\s+=\s+{_NUM}", out)
    return float(m.group(1)) if m else float("nan")


# ----------------------------------------------------------------------
# Figuras
# ----------------------------------------------------------------------
def _png(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def _soma(x, K, EXP):
    s, freq = 0.0, K
    for _ in range(EXP):
        s = s + np.sin(freq * np.pi * x)
        freq *= K
    return s


def _mapa(ax, x, y, Z, titulo, cmap="RdBu_r", log=False):
    if log:
        pos = Z[Z > 0]
        vmax = float(Z.max()) if Z.size else 1.0
        vmin = max(float(pos.min()) if pos.size else 1e-16, vmax * 1e-8, 1e-300)
        m = ax.pcolormesh(x, y, np.ma.masked_less_equal(Z, 0), shading="auto",
                          cmap="inferno", norm=LogNorm(vmin=vmin, vmax=max(vmax, vmin * 10)))
    else:
        m = ax.pcolormesh(x, y, Z, shading="auto", cmap=cmap)
    ax.set_title(titulo, fontsize=10)
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_aspect("equal")
    plt.colorbar(m, ax=ax, fraction=0.046, pad=0.04)


def fig_solucao(pasta, N2, K, EXP, titulo):
    arq = os.path.join(pasta, "grafico.dat")
    if not os.path.isfile(arq) or os.path.getsize(arq) == 0:
        raise RuntimeError("O executável não gerou grafico.dat.")
    try:
        import pandas as pd
        d = pd.read_csv(arq, sep=r"\s+", header=None).to_numpy(dtype=float)
    except ImportError:
        d = np.loadtxt(arq)
    X = d[:, 0].reshape(N2, N2); Y = d[:, 1].reshape(N2, N2)
    U = d[:, 2].reshape(N2, N2); E = d[:, 3].reshape(N2, N2)
    s = max(1, int(np.ceil(N2 / 400)))
    xs, ys = X[::s, 0], Y[0, ::s]
    U, E = U[::s, ::s].T, E[::s, ::s].T
    A = np.outer(_soma(ys, K, EXP), _soma(xs, K, EXP))
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4), constrained_layout=True)
    _mapa(ax[0], xs, ys, A, "Analítica")
    _mapa(ax[1], xs, ys, U, "Numérica")
    _mapa(ax[2], xs, ys, E, "Erro absoluto", log=True)
    fig.suptitle(titulo)
    return _png(fig)


def fig_residuo(runs):
    series = [(r["nx"], [v for v in r["hist"] if math.isfinite(v) and v > 0])
              for r in runs if r["hist"]]
    series = [(n, h) for n, h in series if h]
    if not series:
        return None
    fig, ax = plt.subplots(figsize=(7, 3.4))
    for nx, h in series:
        ax.semilogy(range(1, len(h) + 1), h, "o-", ms=3, label=f"N = {nx}")
    ax.set_xlabel("ciclo V"); ax.set_ylabel("resíduo")
    ax.grid(alpha=.3, which="both"); ax.legend()
    ax.set_title("Histórico do resíduo (após o FMG)", fontsize=10)
    return _png(fig)


def fig_tempos(runs):
    partes = [("chute (C)", "t_guess"), ("pré-suav.", "t_pre"), ("pós-suav.", "t_pos"),
              ("B (grossa)", "t_B"), ("V-cycles (total)", "t_vcycle")]
    fig, ax = plt.subplots(figsize=(8, 3.6))
    larg = 0.8 / max(1, len(runs))
    for j, r in enumerate(runs):
        vals = [1000 * r[k] if not math.isnan(r[k]) else 0 for _, k in partes]
        xs = np.arange(len(partes)) + j * larg
        ax.bar(xs, vals, larg, label=f"N = {r['nx']}")
    ax.set_xticks(np.arange(len(partes)) + larg * (len(runs) - 1) / 2)
    ax.set_xticklabels([n for n, _ in partes], fontsize=9)
    ax.set_ylabel("tempo (ms)"); ax.set_yscale("log"); ax.grid(alpha=.3, axis="y", which="both")
    ax.legend()
    ax.set_title("Onde o tempo é gasto", fontsize=10)
    return _png(fig)


def _fmt(x):
    if isinstance(x, float):
        return "—" if math.isnan(x) else f"{x:.6e}" if (x != 0 and (abs(x) < 1e-3 or abs(x) >= 1e5)) else f"{x:.6g}"
    return str(x)


def tabela_resultados(runs, p_ordem):
    linhas = []
    rot = [("Malha N", "nx"), ("Ciclos V após o FMG", "ciclos_fmg"),
           ("Iterações em C", "iter_C"), ("Resíduo final", "residuo"),
           ("Norma 2 do erro", "erro2"), ("Norma ∞ do erro", "erro_inf"),
           ("Tempo de execução (s)", "tempo"), ("Convergiu", "convergiu")]
    for nome, k in rot:
        v = [("sim" if r[k] else "não") if k == "convergiu" else _fmt(r[k]) for r in runs]
        linhas.append([nome] + v)
    linhas.append(["Ordem de convergência p", _fmt(p_ordem)] + [""] * (len(runs) - 1))
    return linhas


def _corta(texto, cabeca=45, cauda=70):
    L = texto.splitlines()
    if len(L) <= cabeca + cauda + 5:
        return texto
    return "\n".join(L[:cabeca] + [f"... ({len(L) - cabeca - cauda} linhas omitidas) ..."] + L[-cauda:])


# ----------------------------------------------------------------------
# Rotas
# ----------------------------------------------------------------------
@app.route("/")
def index():
    return render_template(
        "index.html", padrao_js=json.dumps(PADRAO), faixas_js=json.dumps(FAIXAS),
        kmax=KMAX, kmax_comparar=KMAX_COMPARAR,
        nomes_js=json.dumps({"A": NOMES_A, "B": NOMES_B, "C": NOMES_C}))


@app.route("/api/executar")
def api_executar():
    try:
        p = _parametros(request.args)
        _, plot = _garante_executaveis(p["EXP"])
        entrada = monta_entrada(p)
        tmp = tempfile.mkdtemp(prefix="lab_")
        try:
            rc, out, err = _roda(plot, p["N"], entrada, tmp)
            runs = parse_execucoes(out)
            if rc != 0 or len(runs) < 2:
                raise RuntimeError((out[-1500:] + err).strip() or "execução falhou")
            ordem = parse_ordem(out)
            if any(not math.isfinite(r["residuo"]) or not math.isfinite(r["erro_inf"])
                   for r in runs):
                raise RuntimeError(
                    "A iteração divergiu (resíduo/erro = NaN ou infinito). Reduza o peso ω "
                    "(Jacobi costuma exigir ω ≤ 1) ou revise ν₁, ν₂ e as tolerâncias.")
            figuras = [("Solução analítica, numérica e erro (malha 2N)",
                        fig_solucao(tmp, 2 * p["N"], p["K"], p["EXP"],
                                    f"Poisson 2D — 2N = {2 * p['N']} — K = {p['K']:g}, EXP = {p['EXP']}"))]
            fr = fig_residuo(runs)
            if fr:
                figuras.append(("Convergência do resíduo", fr))
            figuras.append(("Perfil de tempo", fig_tempos(runs)))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        return jsonify({
            "comando": f"./lab {p['N']} 1",
            "entrada": entrada.replace("\n", " "),
            "resultados": tabela_resultados(runs, ordem),
            "figuras": [{"titulo": t, "img": i} for t, i in figuras],
            "log": _corta(out),
        })
    except ValueError as e:
        return jsonify({"erro": str(e)}), 400
    except subprocess.TimeoutExpired:
        return jsonify({"erro": "Tempo limite excedido (a configuração pode estar divergindo: "
                                "revise os pesos ω, NU1/NU2 e as tolerâncias)."}), 504
    except RuntimeError as e:
        return jsonify({"erro": str(e)}), 500


def _um_combo(exe, p, a, b, c):
    q = dict(p, A=a, B=b, C=c,
             omA=_padrao_omega(a == 2), omB=_padrao_omega(b == 5), omC=_padrao_omega(c == 2))
    tmp = tempfile.mkdtemp(prefix="lab_")
    try:
        rc, out, err = _roda(exe, p["N"], monta_entrada(q), tmp, timeout=60)
    except subprocess.TimeoutExpired:
        return {"a": a, "b": b, "c": c, "erro": "tempo limite"}
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    runs = parse_execucoes(out)
    if rc != 0 or len(runs) < 2:
        return {"a": a, "b": b, "c": c, "erro": "falhou"}
    nv = sum(r["n_vcycle"] for r in runs)
    npre = sum(r["n_pre"] for r in runs)
    npos = sum(r["n_pos"] for r in runs)
    tv = sum(r["t_vcycle"] for r in runs)
    return {
        "a": a, "b": b, "c": c, "p": parse_ordem(out), "vcyc": nv,
        "t_cic_ms": 1000 * tv / nv if nv else 0.0,
        "t_pre_ms": 1000 * sum(r["t_pre"] for r in runs) / npre if npre else 0.0,
        "t_pos_ms": 1000 * sum(r["t_pos"] for r in runs) / npos if npos else 0.0,
        "ncoar": sum(r["n_B_iter"] if b >= 4 else r["n_B_stages"] for r in runs),
        "t": sum(r["tempo"] for r in runs),
    }


def _heat(ax, M, titulo, fmt, log=False, cmap="viridis_r"):
    M = np.array(M, dtype=float)
    ok = M[np.isfinite(M)]
    if log and ok.size and (ok > 0).all():
        im = ax.imshow(M, cmap=cmap, aspect="auto", norm=LogNorm(vmin=ok.min(), vmax=max(ok.max(), ok.min() * 1.0001)))
    else:
        im = ax.imshow(M, cmap=cmap, aspect="auto")
    ax.set_xticks(range(5)); ax.set_xticklabels([NOMES_C[c] for c in range(1, 6)], rotation=30, ha="right", fontsize=8)
    ax.set_yticks(range(5)); ax.set_yticklabels([NOMES_B[b] for b in range(1, 6)], fontsize=8)
    ax.set_xlabel("C (chute inicial)"); ax.set_ylabel("B (malha grossa)")
    ax.set_title(titulo, fontsize=10)
    for i in range(5):
        for j in range(5):
            v = M[i, j]
            ax.text(j, i, "—" if not np.isfinite(v) else fmt(v), ha="center", va="center",
                    fontsize=7.5, color="white" if (np.isfinite(v) and ok.size and
                                                    (v - ok.min()) / max(ok.max() - ok.min(), 1e-300) > 0.5) else "black")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


@app.route("/api/comparar")
def api_comparar():
    """Equivalente ao modo 'esgotar' do lab (50 combinações A/B/C)."""
    try:
        p = _parametros(request.args, comparar=True)
        exe, _ = _garante_executaveis(p["EXP"])
        combos = [(a, b, c) for a in (1, 2) for b in range(1, 6) for c in range(1, 6)]
        with ThreadPoolExecutor(max_workers=4) as ex:
            res = list(ex.map(lambda t: _um_combo(exe, p, *t), combos))
        dados = {(r["a"], r["b"], r["c"]): r for r in res}

        figuras = []
        for chave, titulo, fmt, log in [
            ("t", "Tempo total (s)", lambda v: f"{v:.3g}", True),
            ("vcyc", "Ciclos V (N + 2N)", lambda v: f"{v:.0f}", False),
            ("p", "Ordem de convergência p", lambda v: f"{v:.4f}", False),
            ("ncoar", "Operações na malha grossa", lambda v: f"{v:.0f}", True),
        ]:
            fig, ax = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)
            for k, a in enumerate((1, 2)):
                M = [[dados[(a, b, c)].get(chave, np.nan) for c in range(1, 6)] for b in range(1, 6)]
                _heat(ax[k], M, f"A = {NOMES_A[a]}", fmt, log)
            fig.suptitle(f"{titulo} — N = {p['N']}, K = {p['K']:g}, EXP = {p['EXP']}")
            figuras.append({"titulo": titulo, "img": _png(fig)})

        linhas = []
        for r in sorted(res, key=lambda r: r.get("t", float("inf"))):
            if "erro" in r:
                linhas.append([NOMES_A[r["a"]], NOMES_B[r["b"]], NOMES_C[r["c"]], "erro: " + r["erro"]] + [""] * 6)
            else:
                linhas.append([NOMES_A[r["a"]], NOMES_B[r["b"]], NOMES_C[r["c"]],
                               f"{r['p']:.4f}", str(r["vcyc"]), f"{r['t_cic_ms']:.4f}",
                               f"{r['t_pre_ms']:.4f}", f"{r['t_pos_ms']:.4f}",
                               str(r["ncoar"]), f"{r['t']:.3f}"])
        return jsonify({"cabecalho": ["A", "B", "C", "p", "Vcyc", "t/ciclo (ms)", "tPre (ms)",
                                      "tPos (ms)", "Ncoar", "t (s)"],
                        "linhas": linhas, "figuras": figuras, "N": p["N"]})
    except ValueError as e:
        return jsonify({"erro": str(e)}), 400
    except RuntimeError as e:
        return jsonify({"erro": str(e)}), 500


if __name__ == "__main__":
    if not os.path.isfile(os.path.join(SRC_DIR, "Makefile")):
        raise SystemExit(f"Makefile não encontrado em {SRC_DIR}. Coloque a pasta lab/ "
                         "ao lado de app.py ou defina LAB_DIR.")
    print(f"Fontes: {SRC_DIR}")
    app.run(host="0.0.0.0", port=5000, debug=False, threaded=True)

#!/usr/bin/env python3
"""Figures from results/analysis (only those whose arms were run) -> results/figures/*.png. Exploratory results; zonation profiles are hypotheses."""
import json
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from common import *
cells, vox, seg_h, seg_v = cells_info(); cnt = np.bincount(seg_h, minlength=8); IN = inlet_mM(); A = arms()
D = {a: np.load(ANA / f"ANALYSIS_{a}.npz") for a in A if (ANA / f"ANALYSIS_{a}.npz").exists()}; plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 130})
def seg(a, k, i=-1): return D[a]["segflux"][i][KEY.index(k)]
def have(*names): return all(n in D for n in names)
def save(fig, name): fig.tight_layout(); fig.savefig(FIG / name, dpi=200); plt.close(fig); print("wrote", name)
# ---- prescribed zonation
Z = [("Z0_uniform", "Sem zonação (s=1)", "#555555"), ("Z3_s3_k15", "Zonação s=3", "#1b7837"), ("Z10_s10_k15", "Zonação s=10", "#b2182b")]; Z = [z for z in Z if z[0] in D]
if len(Z) >= 2:
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.4))
    for k, (key, t) in enumerate([("r0123", "PEPCK"), ("r0396", "G6Pase (saída de glicose)"), ("GLYCOGEN_SUPPLY", "Glicogenólise")]):
        for a, lab, c in Z: axs[k].plot(range(8), seg(a, key), "o-", color=c, label=lab)
        axs[k].set_title(t); axs[k].set_xlabel("segmento PP→PC"); axs[k].set_ylabel("fluxo por hospedeiro")
    axs[0].legend(frameon=False); save(fig, "fig_fluxos_por_segmento.png")
    fig, axs = plt.subplots(2, 3, figsize=(12, 6), sharex=True)
    for ax, (f, t) in zip(axs.flat, [("lactate", "Lactato"), ("glycerol", "Glicerol"), ("oxygen", "O₂"), ("ammonia", "Amônia"), ("alanine", "Alanina"), ("palmitate", "Palmitato")]):
        for a, lab, c in Z: v = D[a]["ves"][-1][NAMES.index(f)]; ax.plot(range(8), v / v[0] if f == "glycerol" else v / IN[f], "o-", color=c, label=lab)
        ax.set_title(t + (" (rel. ao segmento 0)" if f == "glycerol" else " (rel. à entrada)")); ax.set_ylim(bottom=0)
    axs[0, 0].legend(frameon=False); save(fig, "fig_perfis_sangue.png")
    mask = np.hypot(*np.meshgrid(-500 + 20 * (np.arange(50) + .5), -500 + 20 * (np.arange(50) + .5))) <= 490
    for f, t, cm in [("lactate", "Lactato (mM)", "viridis"), ("glycerol", "Glicerol (mM)", "viridis"), ("oxygen", "O₂ dissolvido (mM)", "magma")]:
        vals = [np.where(mask, D[a]["maps"][-1][:, NAMES.index(f)].reshape(50, 50), np.nan) for a, _, _ in Z]; lo, hi = np.nanmin(vals), np.nanmax(vals)
        fig, axs = plt.subplots(1, len(Z), figsize=(4.2 * len(Z), 3.8)); axs = np.atleast_1d(axs)
        for ax, (a, l, _), g in zip(axs, Z, vals): im = ax.imshow(g, origin="lower", cmap=cm, vmin=lo, vmax=hi, extent=[-500, 500, -500, 500]); ax.set_title(l); ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(im, ax=axs, shrink=.8, label=t); fig.savefig(FIG / f"fig_mapa_{f}.png", dpi=200, bbox_inches="tight"); plt.close(fig); print("wrote", f"fig_mapa_{f}.png")
# ---- emergent vs prescribed
E = [("Z0_uniform", "Prescrito s=1", "#555555", "--"), ("Z3_s3_k15", "Prescrito s=3", "#1b7837", "--"), ("E0_uniform_k10", "Emergente uniforme κ=1", "#888888", "-"), ("E3_s3_k10", "Emergente s=3 κ=1", "#2166ac", "-"), ("E0_uniform_k15", "Emergente uniforme κ=1,5", "#aaaaaa", "-"), ("E3_s3_k15", "Emergente s=3 κ=1,5", "#b2182b", "-")]; E = [e for e in E if e[0] in D]
if any(e[0].startswith("E") for e in E):
    fig, ax = plt.subplots(1, 3, figsize=(14, 3.8))
    for a, lab, c, ls in E: ax[0].plot(range(8), seg(a, "EX_HC00040_s"), "o" + ls, color=c, label=lab); ax[1].plot(range(8), seg(a, "r0123"), "o" + ls, color=c); ax[2].plot(range(8), D[a]["ves"][-1][NAMES.index("lactate")] / IN["lactate"], "o" + ls, color=c)
    for x, t in zip(ax, ["Exportação de glicose por hospedeiro", "Fluxo de PEPCK por hospedeiro", "Lactato no sangue (rel. à entrada)"]): x.set_title(t); x.set_xlabel("segmento PP→PC")
    ax[0].legend(frameon=False, fontsize=7); save(fig, "fig_emergente_vs_prescrito.png")
# ---- ATP calibration: O2 gradient
if "A2_uniform_pal449" in D:
    def po2(c, s_in=0.95, p50=26.8, n=2.7): S = np.clip(c * s_in, 1e-3, .999); return p50 * (S / (1 - S)) ** (1 / n)
    rel = lambda a: D[a]["ves"][-1][NAMES.index("oxygen")] / IN["oxygen"]; fig, ax = plt.subplots(1, 3, figsize=(14, 3.8))
    for a, lab, c in [("Z0_uniform", "ATPM=5 (antes)", "#888888"), ("A2_uniform_pal449", "ATPM=21,6, enzimas uniformes", "#2166ac"), ("A3_s3_pal449", "ATPM=21,6, zonação s=3", "#b2182b")]:
        if a in D: ax[0].plot(range(8), rel(a), "o-", color=c, label=lab); ax[1].plot(range(8), po2(rel(a)), "o-", color=c); ax[2].plot(range(8), -seg(a, "EX_HC00017_s"), "o-", color=c)
    ax[0].set_title("Conteúdo de O₂ no sangue (rel. à entrada)"); ax[0].legend(frameon=False, fontsize=8); ax[1].axhspan(60, 65, color="#ccc", alpha=.6); ax[1].axhspan(30, 35, color="#ccc", alpha=.6); ax[1].set_title("PO₂ estimada* (mmHg); faixas: lit. 60–65 e 30–35"); ax[2].set_title("Captação de O₂ por hospedeiro")
    for a in ax: a.set_xlabel("segmento PP→PC")
    fig.text(.01, .005, "*curva de Hill externa ao modelo (P50 26,8 mmHg, n 2,7, saturação de entrada 0,95); aproximação.", fontsize=7); fig.tight_layout(rect=(0, .03, 1, 1)); fig.savefig(FIG / "fig_atp_gradiente_O2.png", dpi=200); plt.close(fig); print("wrote fig_atp_gradiente_O2.png")

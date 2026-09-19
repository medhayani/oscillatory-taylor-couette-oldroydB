# -*- coding: utf-8 -*-
"""Apercu de la campagne G pour l'auteur : donnees recues seulement, sans F'.

Lit resultats_G/G_S<S>_E<E>.csv (consolides par _G_quicklook.py).
Une figure : colonnes E = 0.01, 1, 10, 50 ; lignes Wi_c(gamma), k_c(gamma) ;
une courbe par S.
Zones grisees : De = 2 E gamma^2 < 0.3, ou l'operateur de la campagne G n'a
pas le terme de Coriolis azimutal (annexe A.5, _epsv_batch.log) :
  gris clair 0.1 <= De < 0.3 (ecart 2 a 20 %), gris fonce De < 0.1 (>= 20 %,
  jusqu'a x69). Trait pale dans ces zones. Le minimum est un disque plein s'il
  tombe hors zone grisee, un cercle vide sinon.
Palette (dataviz/validate_palette.js, clair) : ALL PASS, contraste WARN pour
aqua et jaune -> etiquettes directes S = ... en bout de courbe.
Sortie : resultats_G/G_apercu.png ; minima : resultats_G/G_apercu_minima.txt.
Usage : python _G_apercu.py
"""
import os, json, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import LogLocator, FuncFormatter, NullFormatter

RG = "resultats_G"
SQ = np.sqrt(0.14)
SS = (0.98, 0.9, 0.7, 0.3)
ES = (0.01, 1.0, 10.0, 50.0)
COL = {0.98: "#2a78d6", 0.9: "#eb6834", 0.7: "#1baf7a", 0.3: "#eda100"}
SURF, INK1, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BAND1, BAND2 = "#efeeea", "#dcdbd4"      # 0.1 <= De < 0.3 ; De < 0.1
GAP = 1.35
X0, X1 = 0.009, 22.0


def recu():
    """Date de la derniere donnee consolidee (et non une heure ecrite a la main)."""
    fs = [os.path.join(RG, f) for f in os.listdir(RG) if f.startswith("G_S") and f.endswith(".csv")]
    return time.strftime("%d/%m à %H:%M", time.localtime(max(os.path.getmtime(f) for f in fs)))


def couverture():
    st = json.load(open("G_state.json"))
    ids = json.load(open(os.path.join("kg_G", "ids.json")))
    out = []
    for pre, nom in (("shqgl-", "basse fréquence"), ("shqg-", "corps"), ("shqgt-", "queues")):
        ks = [k for k in ids if k.startswith(pre)]
        d = sum(1 for k in ks if st.get(k, {}).get("done"))
        out.append("%s %d/%d" % (nom, d, len(ks)))
    return "Kernels terminés : " + ", ".join(out) + "."


def load(S, E):
    f = os.path.join(RG, "G_S%.2f_E%g.csv" % (S, E))
    if not os.path.exists(f):
        return None
    d = np.atleast_2d(np.genfromtxt(f, delimiter=",", skip_header=1))
    d = d[np.isfinite(d[:, 1])]
    return (d[:, 0], d[:, 3], d[:, 2]) if len(d) else None


def gaps(g, *ys):
    cut = np.flatnonzero(g[1:] / g[:-1] > GAP) + 1
    return (np.insert(g, cut, np.nan),) + tuple(np.insert(y, cut, np.nan) for y in ys)


def main():
    plt.rcParams.update({
        "font.size": 8.5, "font.family": "DejaVu Sans",
        "axes.edgecolor": AXIS, "axes.linewidth": 0.7, "axes.labelcolor": INK2,
        "xtick.color": INK2, "ytick.color": INK2, "text.color": INK1,
        "lines.solid_capstyle": "round", "lines.solid_joinstyle": "round"})
    fig, axes = plt.subplots(2, 4, figsize=(14.0, 7.0), sharex=True, facecolor=SURF,
                             gridspec_kw={"height_ratios": [1.45, 1.0], "hspace": 0.1, "wspace": 0.4})
    rows = []
    for j, E in enumerate(ES):
        a1, a2 = axes[:, j]
        g01, g03 = np.sqrt(0.1 / (2 * E)), np.sqrt(0.3 / (2 * E))
        for ax in (a1, a2):
            ax.set_facecolor(SURF)
            ax.set_xscale("log")
            ax.set_xlim(X0, X1)
            if g03 > X0:
                ax.axvspan(X0, min(g03, X1), color=BAND1, lw=0, zorder=0)
            if g01 > X0:
                ax.axvspan(X0, min(g01, X1), color=BAND2, lw=0, zorder=0)
            ax.grid(True, which="major", color=GRID, lw=0.6, zorder=1)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
        a1.set_yscale("log")
        a1.yaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
        a1.yaxis.set_major_formatter(FuncFormatter(lambda y, _: "%g" % y))
        a1.yaxis.set_minor_formatter(NullFormatter())
        a2.xaxis.set_major_formatter(FuncFormatter(lambda x, _: "%g" % x))
        missing, ends = [], []
        for S in SS:
            d = load(S, E)
            if d is None:
                missing.append(S)
                continue
            g, wi, k = d
            c = COL[S]
            ok = 2 * E * g ** 2 >= 0.3
            for sel, alpha, lw in ((~ok, 0.35, 1.3), (ok, 1.0, 1.8)):
                if sel.sum() >= 2:
                    gg, ww, kk = gaps(g[sel], wi[sel], k[sel])
                    a1.plot(gg, ww, color=c, lw=lw, alpha=alpha, zorder=3)
                    a2.plot(gg, kk, color=c, lw=lw, alpha=alpha, zorder=3)
            i = int(np.argmin(wi))
            fiable = 2 * E * g[i] ** 2 >= 0.3
            a1.plot([g[i]], [wi[i]], "o", ms=6, mfc=c if fiable else SURF, mec=c, mew=1.6, zorder=5)
            ends.append([np.log10(wi[-1]), g[-1], "S=%g" % S])
            rows.append((E, S, wi[i], g[i], k[i], 2 * E * g[i] ** 2, fiable, g.max()))
        # etiquettes directes en bout de courbe, ecartees d'au moins 5 % de la hauteur (echelle log)
        lo, hi = (np.log10(v) for v in a1.get_ylim())
        sep = 0.055 * (hi - lo)
        ends.sort()
        for n in range(1, len(ends)):
            ends[n][0] = max(ends[n][0], ends[n - 1][0] + sep)
        for ly, gx, lab in ends:
            a1.text(X1 * 1.05, 10 ** ly, lab, color=INK2, fontsize=7.5, va="center", ha="left", clip_on=False)
        title = "E = %g" % E + ("   (S = %s : en cours)" % ", ".join("%g" % s for s in missing) if missing else "")
        a1.set_title(title, color=INK1, fontsize=11, pad=8)
        a2.set_xlabel("γ", fontsize=10)
    axes[0, 0].set_ylabel("Wi_c", fontsize=10)
    axes[1, 0].set_ylabel("k_c", fontsize=10)
    fig.suptitle("Campagne G (inertie complète) : données reçues le %s" % recu(),
                 x=0.055, y=0.985, ha="left", fontsize=12.5, color=INK1, fontweight="bold")
    fig.text(0.055, 0.935, couverture() + " "
             "Zones grisées : terme de Coriolis manquant, seuil non fiable. "
             "Disque plein : minimum hors zone grisée ; cercle vide : minimum dans la zone grisée.",
             ha="left", fontsize=8.5, color=INK2)
    handles = [Line2D([0], [0], color=COL[S], lw=2.2, label="S = %g%s" % (S, "  (pour S = 1)" if S == 0.98 else ""))
               for S in SS]
    handles += [Patch(facecolor=BAND1, edgecolor=AXIS, lw=0.5, label="0.1 ≤ De < 0.3 : écart 2 à 20 %"),
                Patch(facecolor=BAND2, edgecolor=AXIS, lw=0.5, label="De < 0.1 : écart ≥ 20 % (jusqu'à ×69)")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=6,
               frameon=False, fontsize=9)
    fig.subplots_adjust(left=0.055, right=0.965, top=0.87, bottom=0.14)
    out = os.path.join(RG, "G_apercu.png")
    fig.savefig(out, dpi=140, facecolor=SURF)
    plt.close(fig)
    lines = ["E      S      min Wi_c   gamma*   k_c   De*      fiable   gamma max recu"]
    for E, S, w, gs, kc, de, fi, gm in rows:
        lines.append("%-6g %-6g %8.4g  %7.3f  %4g  %-7.3g  %-7s  %.2f" % (E, S, w, gs, kc, de, "oui" if fi else "non", gm))
    txt = "\n".join(lines)
    open(os.path.join(RG, "G_apercu_minima.txt"), "w", encoding="utf-8").write(txt + "\n")
    print(txt)
    print("figure :", out)


if __name__ == "__main__":
    main()

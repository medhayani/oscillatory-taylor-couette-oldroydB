# -*- coding: utf-8 -*-
"""Figures de la campagne G : UNE figure par cas (S, E), sans titre.

Rendu LaTeX reel (newtxtext/newtxmath, comme les figures de la campagne F'),
axes et etiquettes en grand. Wi_c(gamma) en noir (log-log, axe de gauche),
k_c(gamma) en bleu (axe de droite, lineaire). Pas de marque au minimum
(demande de l'auteur, 19/09) : les valeurs sont dans minima.txt.
k_c est trace en tirets la ou il touche le bord de la grille en k (40 sous
gamma = 7.8, 80 au-dessus) : le minimum en k n'y est pas atteint et Wi_c y est
un majorant. Pas de bande grise, pas de loi WKB en gamma^-2 (annexe C.3 : elle
est la limite du systeme sans inertie de la perturbation, pas de (P1)-(P8)).

Sorties : resultats_G/figures_G/G_E<tag>_S<tag>.png (16 figures, 300 dpi),
copiees dans article_JFM/figures/ (figures 2 a 5 de l'article), et
resultats_G/figures_G/minima.txt.
Si LaTeX n'est pas disponible, bascule automatiquement sur mathtext.
Usage : python _G_figs_article.py
"""
import os
import shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RG = "resultats_G"
OUT = os.path.join(RG, "figures_G")
FIG_ART = os.path.join("article_JFM", "figures")
SS = (0.98, 0.9, 0.7, 0.3)
ES = (0.01, 1.0, 10.0, 50.0)
C_WI, C_K, GRID, C_WKB = "#000000", "#1f5fa9", "#e3e2dc", "#c0392b"
# Loi sans inertie (annexe C.3) : Wi_c = eps^{-1/2} A(S) / gamma^2, avec A(S)
# de article_JFM/wkb_table.tex (campagne quasi statique). Tracee en tirets
# rouges sur la branche WKB, De = 2 E gamma^2 <= 0.05, comme repere : le
# systeme complet a E fixe ne la suit pas (pente mesuree -1, non -2).
A_WKB = {0.98: 5.50, 0.9: 5.88, 0.7: 6.99, 0.3: 11.47}
X0, X1 = 0.009, 22.0
GAP = 1.35
LOG = True           # True : axes logarithmiques ; False : axes lineaires
KMAX_BAS, KMAX_QUEUE = 40.0, 80.0

STYLE_LATEX = {
    "text.usetex": True,
    "font.family": "serif",
    "text.latex.preamble": r"\usepackage{newtxtext,newtxmath}",
    "font.size": 16, "axes.labelsize": 21,
    "xtick.labelsize": 16, "ytick.labelsize": 16,
    "axes.linewidth": 1.0, "axes.edgecolor": "#222222",
    "xtick.major.width": 1.0, "ytick.major.width": 1.0,
    "xtick.major.size": 5, "ytick.major.size": 5,
}
STYLE_SANS_LATEX = dict(STYLE_LATEX, **{"text.usetex": False,
                                        "mathtext.fontset": "stix"})


def tag(v):
    return ("%g" % v).replace(".", "p")


def load(S, E):
    f = os.path.join(RG, "G_S%.2f_E%g.csv" % (S, E))
    if not os.path.exists(f):
        return None
    d = np.atleast_2d(np.genfromtxt(f, delimiter=",", skip_header=1))
    d = d[np.isfinite(d[:, 1])]
    return (d[:, 0], d[:, 3], d[:, 2]) if len(d) else None


def complet(d):
    """Cas calcule a 100 % : 19 points de basse frequence, 375 de corps, 62 de queue."""
    g = d[0]
    return ((g < 0.1).sum() >= 19 and ((g >= 0.1) & (g < 7.78)).sum() >= 375
            and (g >= 7.78).sum() >= 62)


def kmax_of(g):
    """Bord de la grille en k, selon la famille dont vient le point."""
    return np.where(g < 7.78, KMAX_BAS, KMAX_QUEUE)


def cut(g, *ys):
    """Coupe le trait la ou la couverture en gamma a un trou."""
    if len(g) < 2:
        return (g,) + ys
    c = np.flatnonzero(g[1:] / g[:-1] > GAP) + 1
    return (np.insert(g, c, np.nan),) + tuple(np.insert(y, c, np.nan) for y in ys)


def dessine(d, chemin, S=None, E=None):
    g, wi, k = d
    fig, ax = plt.subplots(figsize=(6.4, 4.6), facecolor="white")
    if LOG:
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(X0, X1)
        ax.set_ylim(float(np.nanmin(wi)) / 1.6, float(np.nanmax(wi)) * 1.6)   # cadre sur les donnees
    else:                                   # axes lineaires (demande du 17/09)
        ax.set_xlim(0.0, 20.5)
        ax.set_ylim(0.0, float(np.nanmax(wi)) * 1.05)
    ax.grid(True, which="major", color=GRID, lw=0.6)
    ax.set_xlabel(r"$\gamma$")
    ax.set_ylabel(r"$Wi_c$")
    axk = ax.twinx()
    axk.set_ylabel(r"$k_c$", color=C_K)
    axk.tick_params(axis="y", colors=C_K)
    axk.set_ylim(0, 85)

    if S in A_WKB and E is not None:            # loi sans inertie, branche WKB
        g_hi = min(np.sqrt(0.025 / E), X1)
        if g_hi > X0 * 1.5:
            gw = np.geomspace(X0 * 1.1, g_hi, 40)
            ax.plot(gw, A_WKB[S] / np.sqrt(0.14) / gw ** 2, ls=(0, (5, 3)),
                    lw=2.4, color=C_WKB, zorder=2)
    gg, ww, kk = cut(g, wi, k)
    ax.plot(gg, ww, "-", color=C_WI, lw=1.8, zorder=4)
    axk.plot(gg, kk, "-", color=C_K, lw=1.2, zorder=3)
    sat = k >= kmax_of(g) - 1e-9
    if sat.any():
        gs, kss = cut(g, np.where(sat, k, np.nan))
        axk.plot(gs, kss, ls=(0, (4, 2)), color=C_K, lw=2.0, zorder=5)
    i = int(np.argmin(wi))          # garde pour le tableau des minima, plus marque sur la figure
    fig.tight_layout()
    fig.savefig(chemin, dpi=300, facecolor="white")
    plt.close(fig)
    return i


def main():
    os.makedirs(OUT, exist_ok=True)
    plt.rcParams.update(STYLE_LATEX)
    lines = ["E      S      min Wi_c   gamma*   k_c   De*      queues/62"]
    premier = True
    for E in ES:
        for S in SS:
            d = load(S, E)
            if d is None:
                print("E=%g S=%g : pas encore de donnees" % (E, S))
                continue
            nom = "G_E%s_S%s.png" % (tag(E), tag(S))
            if not complet(d):        # cas non calcule a 100 % : pas de figure
                for f in (os.path.join(OUT, nom), os.path.join(FIG_ART, nom)):
                    if os.path.exists(f):
                        os.remove(f)
                g = d[0]
                print("E=%-5g S=%-5g : incomplet (bf %d/19, corps %d/375, queues %d/62),"
                      " figure retiree" % (E, S, (g < 0.1).sum(),
                                           ((g >= 0.1) & (g < 7.78)).sum(), (g >= 7.78).sum()))
                continue
            f = os.path.join(OUT, nom)
            try:
                i = dessine(d, f, S, E)
            except RuntimeError as e:          # LaTeX absent ou en echec
                if not premier:
                    raise
                print("LaTeX indisponible (%s) -> rendu mathtext" % str(e).splitlines()[0][:60])
                plt.rcParams.update(STYLE_SANS_LATEX)
                i = dessine(d, f, S, E)
            premier = False
            shutil.copy(f, os.path.join(FIG_ART, nom))
            g, wi, k = d
            lines.append("%-6g %-6g %8.4g  %7.3f  %4g  %-7.3g  %d"
                         % (E, S, wi[i], g[i], k[i], 2 * E * g[i] ** 2, int((g >= 7.78).sum())))
            print("figure :", f)
    txt = "\n".join(lines)
    open(os.path.join(OUT, "minima.txt"), "w", encoding="utf-8").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()

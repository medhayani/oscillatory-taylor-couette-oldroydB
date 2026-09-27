# -*- coding: utf-8 -*-
"""Tableau 3 de l'article (tab:shq) a partir des courbes assemblees GH.

Pour chaque cas : gamma* et De* du minimum de Wi_c, Wi_c* et k_c*, plus les
deux seuils aux bords, Wi_c(0.1) et Wi_c(20), et le rapport entre eux, qui
dit ce que la modulation gagne. Ecrit article_JFM/shq_table.tex (copie par
_build_prf.py). Un cas dont la courbe est incomplete ou en dents de scie sur
[0.1, 20] a ses valeurs quand meme : le minimum n'est jamais dans la zone
concernee (gamma <= 0.4 a E = 0.01, queues de S = 1).

Usage : python _GH_table.py
"""
import os

import numpy as np

RG = "resultats_G"
SS = (0.98, 0.9, 0.7, 0.3)
ES = (0.01, 1.0, 10.0, 50.0)
LAB = {0.98: "1", 0.9: "0.9", 0.7: "0.7", 0.3: "0.3"}


def fmt(v, n=3):
    if v >= 1000:
        m, e = v / 10 ** int(np.floor(np.log10(v))), int(np.floor(np.log10(v)))
        return r"$%.2f\times10^{%d}$" % (m, e)
    return "$%s$" % ("%.*g" % (n, v))


def main():
    lignes = ["% genere par _GH_table.py (courbes assemblees GH, 23/09/2026)",
              r"\begin{tabular}{cccccccc}", r"\toprule",
              r"$S$ & $E$ & $\gamma^*$ & $\De^{*}$ & $\Wi_c^{*}$ & $k_c^{*}$ & "
              r"$\Wi_c(0.2)$ & $\Wi_c(20)$\\", r"\midrule"]
    print("S     E      gamma*  De*     Wi_c*   k_c*   Wi(0.2)  Wi(20)")
    for E in ES:
        for S in SS:
            d = np.genfromtxt(os.path.join(RG, "GH_S%.2f_E%g.csv" % (S, E)), delimiter=",", skip_header=1)
            # Meme borne basse que les figures : a E = 0.01 les points sous
            # gamma = 0.18 ne sont ni dans le domaine de l'entrefer etroit ni
            # converges a N = 28 (voir G_MIN dans _G_figs_article.py).
            d = d[np.isfinite(d[:, 1]) & (d[:, 0] >= (0.18 if E < 0.1 else 0.1) - 1e-9)]
            g, wi, k = d[:, 0], d[:, 3], d[:, 2]
            i = int(np.argmin(wi))
            w01 = wi[int(np.argmin(np.abs(g - 0.2)))]
            w20 = wi[int(np.argmin(np.abs(g - 20.0)))]
            de = 2 * E * g[i] ** 2
            lignes.append("$%s$ & $%g$ & %s & %s & %s & $%g$ & %s & %s\\\\" % (
                LAB[S], E, fmt(g[i]), fmt(de, 2), fmt(wi[i]), k[i], fmt(w01), fmt(w20)))
            print("%-5g %-6g %-7.3g %-7.2g %-7.3g %-6g %-8.3g %-8.3g" % (S, E, g[i], de, wi[i], k[i], w01, w20))
    lignes += [r"\bottomrule", r"\end{tabular}"]
    out = os.path.join("article_JFM", "shq_table.tex")
    open(out, "w", encoding="utf-8").write("\n".join(lignes) + "\n")
    print("ecrit", out)


if __name__ == "__main__":
    main()

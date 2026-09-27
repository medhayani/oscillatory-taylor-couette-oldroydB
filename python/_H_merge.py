# -*- coding: utf-8 -*-
"""Campagne H : fusion des blocs corriges et comparaison a la campagne G.

Lit resultats_G/corr_S<S>_E<E>_p*.csv (paire de Coriolis complete) et ecrit
resultats_G/H_S<S>_E<E>.csv (gamma, Wp_c, k_c, Wi_c). Affiche, pour chaque
cas : la couverture, le rapport H/G au meme gamma, et le produit Wp_c gamma^2
aux trois plus basses frequences, a comparer a A(S) du tableau sans inertie
(5.50, 5.88, 6.99, 11.47 pour S = 1, 0.9, 0.7, 0.3).

Usage : python _H_merge.py
"""
import glob
import os

import numpy as np

RG = "resultats_G"
SQ = np.sqrt(0.14)
SS = (0.98, 0.9, 0.7, 0.3)
ES = (0.01, 1.0, 10.0, 50.0)
A_WKB = {0.98: 5.50, 0.9: 5.88, 0.7: 6.99, 0.3: 11.47}


def lire_blocs(S, E):
    rows = {}
    for f in sorted(glob.glob(os.path.join(RG, "corr_S%.2f_E%g_p*.csv" % (S, E)))):
        for line in open(f).read().splitlines()[1:]:
            v = [x.strip() for x in line.split(",")]
            if len(v) < 3 or not v[0]:
                continue
            w = float(v[1]) if v[1] != "nan" else np.nan
            rows[round(float(v[0]), 5)] = (w, float(v[2]) if v[2] != "nan" else np.nan)
    return rows


def charge_G(S, E):
    f = os.path.join(RG, "G_S%.2f_E%g.csv" % (S, E))
    d = np.genfromtxt(f, delimiter=",", skip_header=1)
    d = d[np.isfinite(d[:, 1])]
    return d[:, 0], d[:, 1]


def main():
    print("cas        points   gamma          H/G au plus bas   Wp_c*gamma^2 (3 plus bas)   A(S)")
    for S in SS:
        for E in ES:
            rows = lire_blocs(S, E)
            if not rows:
                continue
            g = np.array(sorted(rows))
            w = np.array([rows[x][0] for x in g])
            k = np.array([rows[x][1] for x in g])
            out = os.path.join(RG, "H_S%.2f_E%g.csv" % (S, E))
            with open(out, "w") as fh:
                fh.write("gamma,Wp_c,k_c,Wi_c\n")
                for a, b, c in zip(g, w, k):
                    fh.write("%.5f,%.6g,%g,%.6g\n" % (a, b, c, b / SQ))
            ok = np.isfinite(w)
            gG, wG = charge_G(S, E)
            rap = np.nan
            if ok.any():
                i = int(np.flatnonzero(ok)[0])
                j = int(np.argmin(np.abs(gG - g[i])))
                rap = w[i] / wG[j]
            prod = [w[n] * g[n] ** 2 for n in np.flatnonzero(ok)[:3]]
            print("S=%-5g E=%-5g %4d  %.3f-%-6.2f  x%-6.2f  %s   %.2f"
                  % (S, E, int(ok.sum()), g[0], g[-1], rap,
                     "  ".join("%.2f" % p for p in prod), A_WKB[S]))
    print("\nfichiers ecrits : resultats_G/H_S<S>_E<E>.csv")


if __name__ == "__main__":
    main()

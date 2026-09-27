# -*- coding: utf-8 -*-
"""Courbes finales : campagne H (corrigee) sous De = 0.3, campagne G au-dessus.

La campagne G ecarte la moitie azimutale du couplage de Coriolis ; l'ecart au
seuil corrige est sous 2-3 % pour De >= 0.3 et explose en dessous. La campagne
H a recalcule De < 0.3 avec le terme restitue. Ce script assemble, pour chaque
cas, resultats_G/GH_S<S>_E<E>.csv :
    gamma <= dernier gamma de H  -> points de H (les NaN de H sont omis)
    gamma >  dernier gamma de H  -> points de G
et affiche le saut au raccord, a surveiller (attendu : quelques %).

Usage : python _GH_merge.py
"""
import os

import numpy as np

RG = "resultats_G"
SQ = np.sqrt(0.14)
SS = (0.98, 0.9, 0.7, 0.3)
ES = (0.01, 1.0, 10.0, 50.0)


def lit(nom):
    d = np.atleast_2d(np.genfromtxt(os.path.join(RG, nom), delimiter=",", skip_header=1))
    return d[np.isfinite(d[:, 1])]


def main():
    print("cas            H (pts, plage)        raccord   saut H/G   total")
    for S in SS:
        for E in ES:
            g = lit("G_S%.2f_E%g.csv" % (S, E))
            fh_ = os.path.join(RG, "H_S%.2f_E%g.csv" % (S, E))
            if os.path.exists(fh_):
                h = lit("H_S%.2f_E%g.csv" % (S, E))
            else:
                h = np.zeros((0, 4))
            if len(h):
                gmax = h[-1, 0]
                haut = g[g[:, 0] > gmax + 1e-9]
                tout = np.vstack([h, haut])
                j = int(np.argmin(np.abs(g[:, 0] - gmax)))
                saut = h[-1, 1] / g[j, 1]
                info = "%3d, %.3f-%.2f" % (len(h), h[0, 0], gmax)
            else:
                tout, saut, gmax, info = g, np.nan, np.nan, "aucun"
            with open(os.path.join(RG, "GH_S%.2f_E%g.csv" % (S, E)), "w") as fh:
                fh.write("gamma,Wp_c,k_c,Wi_c\n")
                for a, b, c, _ in tout:
                    fh.write("%.5f,%.6g,%g,%.6g\n" % (a, b, c, b / SQ))
            print("S=%-5g E=%-5g %-20s  g=%-5.2f   x%.3f     %d" % (S, E, info, gmax, saut, len(tout)))


if __name__ == "__main__":
    main()

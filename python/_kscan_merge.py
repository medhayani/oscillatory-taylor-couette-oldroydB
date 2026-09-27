# -*- coding: utf-8 -*-
"""Fusion du balayage complet en k (kernels shqk-*) dans les courbes assemblees.

Lit resultats_G/kscan_S<S>_E<E>_p*.csv (gamma, k, Wp_c), prend a chaque gamma
le minimum sur k, et REMPLACE le point correspondant de GH_S<S>_E<E>.csv
(sauvegarde GH_..._avant_kscan.csv la premiere fois). Un gamma sans seuil fini
est laisse tel quel. Affiche, par cas, le nombre de points remplaces et le
plus grand changement.

Usage : python _kscan_merge.py
"""
import glob
import os
import shutil

import numpy as np

RG = "resultats_G"
SQ = np.sqrt(0.14)


def main():
    for f in sorted(glob.glob(os.path.join(RG, "kscan_S*_E*_p*.csv"))):
        pass
    cas = {}
    for f in glob.glob(os.path.join(RG, "kscan_S*_E*_p*.csv")):
        nom = os.path.basename(f)
        S = float(nom.split("_S")[1].split("_E")[0])
        E = float(nom.split("_E")[1].split("_p")[0])
        d = np.atleast_2d(np.genfromtxt(f, delimiter=",", skip_header=1))
        if d.size == 0:
            continue
        cas.setdefault((S, E), []).append(d)
    print("cas            gammas   remplaces   pire changement (gamma, avant -> apres, k)")
    for (S, E), blocs in sorted(cas.items()):
        dall = np.vstack(blocs)          # avec les NaN : sert a compter les k faits
        d = dall[np.isfinite(dall[:, 2])]
        gh_f = os.path.join(RG, "GH_S%.2f_E%g.csv" % (S, E))
        bak = os.path.join(RG, "GH_S%.2f_E%g_avant_kscan.csv" % (S, E))
        if not os.path.exists(bak):
            shutil.copy2(gh_f, bak)
        gh = np.genfromtxt(bak, delimiter=",", skip_header=1)
        pire, nrep, npart, nbord = (None, 1.0), 0, 0, 0
        # nombre de k attendus : famille 1 (E = 0.01) balaie k = 1..80 pas 1,
        # famille 2 (queues S = 0.98) k = 10..80 pas 2. Un gamma coupe en cours
        # de route a une liste tronquee : son "minimum" serait faux, on le saute.
        nk_att = 80 if E < 0.1 else 36
        for gam in np.unique(np.round(dall[:, 0], 5)):
            ma = np.abs(dall[:, 0] - gam) < 1e-6
            if np.unique(np.round(dall[ma, 1], 3)).size < 0.95 * nk_att:
                npart += 1
                continue
            m = np.abs(d[:, 0] - gam) < 1e-6
            if not m.any():
                continue
            i = int(np.argmin(d[m, 2]))
            k_best, w_best = d[m, 1][i], d[m, 2][i]
            # Le minimum doit etre INTERIEUR a la grille en k. S'il tombe au
            # bord, le vrai minimum est au-dela et la valeur n'est pas un
            # seuil : on garde le point de la campagne (23/09).
            if k_best >= np.max(d[m, 1]) - 1e-9 or k_best <= np.min(d[m, 1]) + 1e-9:
                nbord += 1
                continue
            j = int(np.argmin(np.abs(gh[:, 0] - gam)))
            if abs(gh[j, 0] - gam) > 1e-4:
                continue
            fac = w_best / gh[j, 1] if np.isfinite(gh[j, 1]) and gh[j, 1] > 0 else np.nan
            if np.isfinite(fac) and abs(np.log(fac)) > abs(np.log(pire[1])):
                pire = ((gam, gh[j, 1], w_best, k_best), fac)
            gh[j, 1], gh[j, 2], gh[j, 3] = w_best, k_best, w_best / SQ
            nrep += 1
        with open(gh_f, "w") as fh:
            fh.write("gamma,Wp_c,k_c,Wi_c\n")
            for a, b, c, e in gh:
                fh.write("%.5f,%.6g,%g,%.6g\n" % (a, b, c, e))
        p = pire[0]
        print("S=%-5g E=%-5g %4d     %4d (%d partiels, %d au bord)  %s" % (
            S, E, len(np.unique(np.round(dall[:, 0], 5))), nrep, npart, nbord,
            ("gamma=%.2f : %.4g -> %.4g (k=%g), x%.3f" % (p[0], p[1], p[2], p[3], pire[1])) if p else "-"))


if __name__ == "__main__":
    main()

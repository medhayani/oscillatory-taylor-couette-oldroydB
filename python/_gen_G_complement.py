# -*- coding: utf-8 -*-
"""Kernels de rattrapage : fréquences de queue manquantes de la campagne G.

Les queues de S = 0.98 (E = 1, 10 et 50) ont été coupées par la limite de 12 h
de Kaggle : 36, 26 et 24 points sur 62. Les journaux ne contiennent rien de
plus que les CSV (vérifié par _G_recup_logs.py), donc il faut recalculer les
points manquants.

Ce script reprend le corps de _gen_G_kaggle.py (même opérateur, N = 48,
k dans [2, 80]) et n'écrit que les gamma absents de resultats_G/G_S<S>_E<E>.csv,
en kernels courts : NSUB = 4 sous-blocs parallèles, au plus GMAX gamma par
sous-bloc, soit environ 8 h par kernel au rythme mesuré (1.5 à 2 h par point).
Amorçage : (k_c, Wp_c) du point de la campagne G le plus proche, et non de la
campagne F'.

Les kernels sont écrits dans kg_G/<id>/ et leurs identifiants ajoutés à la fin
de kg_G/ids.json : le pilote _G_autopilot.py les pousse à son cycle suivant.
Sorties : tail_S<S>_E<E>_p21.csv et suivants, fusionnés par _G_quicklook.py.

Usage : python _gen_G_complement.py [--ecrire]
        (sans --ecrire : n'affiche que ce qui serait fait)
"""
import json
import os
import sys

import numpy as np

import _gen_G_kaggle as G

RG = "resultats_G"
CAS = ((0.98, 1.0), (0.98, 10.0), (0.98, 50.0))
GMAX = 4                      # gamma par sous-bloc (4 sous-blocs par kernel)
BLOC0 = 21                    # numéros de blocs de sortie, après p1..p8
PREFIXE = "shqtc"             # tail complement
# Generation : les kernels de la 1re vague gardent leur nom ; --gen 2 en cree
# de nouveaux (shqtc2-...) pour les gamma encore manquants, avec d'autres
# numeros de blocs de sortie.
GEN = 1


def manquants(S, E):
    f = os.path.join(RG, "G_S%.2f_E%g.csv" % (S, E))
    d = np.atleast_2d(np.genfromtxt(f, delimiter=",", skip_header=1))
    vus = set(np.round(d[np.isfinite(d[:, 1])][:, 0], 2))
    return [float(g) for g in G.TAIL_G if round(float(g), 2) not in vus]


def amorce(S, E, gam):
    """(k_c, Wp_c) du point de la campagne G le plus proche de gam."""
    d = np.atleast_2d(np.genfromtxt(os.path.join(RG, "G_S%.2f_E%g.csv" % (S, E)),
                                    delimiter=",", skip_header=1))
    d = d[np.isfinite(d[:, 1])]
    i = int(np.argmin(np.abs(d[:, 0] - gam)))
    return float(d[i, 2]), float(d[i, 1])


def main():
    global GEN, GMAX, BLOC0, CAS
    ecrire = "--ecrire" in sys.argv
    if "--gen" in sys.argv:
        GEN = int(sys.argv[sys.argv.index("--gen") + 1])
        BLOC0 = 21 + 20 * (GEN - 1)
    if "--gmax" in sys.argv:
        GMAX = int(sys.argv[sys.argv.index("--gmax") + 1])
    if "--cas" in sys.argv:                 # ex. --cas 0.98,50
        S, E = (float(x) for x in sys.argv[sys.argv.index("--cas") + 1].split(","))
        CAS = ((S, E),)
    corps = G.build_body(ncheb=G.NCHEB_T, a_hi=G.K_HI_T)
    ids_path = os.path.join(G.ROOT, "ids.json")
    ids = json.load(open(ids_path))
    neufs = []
    for S, E in CAS:
        miss = manquants(S, E)
        if not miss:
            print("S=%g E=%g : rien à compléter" % (S, E))
            continue
        par_kernel = G.NSUB * GMAX
        nk = int(np.ceil(len(miss) / par_kernel))
        print("S=%g E=%g : %d gamma manquants -> %d kernels" % (S, E, len(miss), nk))
        for c, part in enumerate(np.array_split(np.array(miss), nk)):
            subs = [list(map(float, x)) for x in np.array_split(part, G.NSUB) if len(x)]
            kid = "%s%s-s%se%s-c%d" % (PREFIXE, "" if GEN == 1 else GEN, G.tag(S), G.tag(E), c + 1)
            seed = {"%.5f" % s[0]: amorce(S, E, s[0]) for s in subs}
            src = corps + G.FOOTER_LIST % (
                repr(S), repr(E), repr(subs), BLOC0 + c * G.NSUB,
                "tail_S%.2f_E%g_p%d.csv", repr(seed))
            print("   %-22s %2d gamma : %s" % (kid, sum(len(s) for s in subs),
                                               ", ".join("%.1f" % x for s in subs for x in s)))
            if ecrire:
                G.write_kernel(kid, src)
            neufs.append(kid)
    if not ecrire:
        print("\nessai a blanc ; relancer avec --ecrire pour ecrire les kernels")
        return
    ids += [k for k in neufs if k not in ids]
    tmp = ids_path + ".tmp"
    json.dump(ids, open(tmp, "w"), indent=1)
    os.replace(tmp, ids_path)                  # remplacement atomique : le pilote lit ids.json
    print("\n%d kernels ecrits, ids.json passe a %d entrees" % (len(neufs), len(ids)))


if __name__ == "__main__":
    main()

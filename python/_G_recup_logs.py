# -*- coding: utf-8 -*-
"""Recupere les points des kernels coupes par Kaggle a 12 h.

run_full_sweep.march imprime chaque point resolu ("ROW S=... E=... gamma, Wp_c,
k_c, Wi_c"). Les journaux rapatries dans kg_G_out/<kernel>/*.log sont du JSON
Kaggle : chaque entree porte ce texte dans son champ "data". On en extrait les
points et on ecrit un bloc supplementaire resultats_G/tail_S<S>_E<E>_plog.csv,
que _G_quicklook.py fusionne comme les autres blocs.

Usage : python _G_recup_logs.py
"""
import glob
import json
import os
import re

RG = "resultats_G"
KO = "kg_G_out"
ROW = re.compile(r"ROW S=(\S+) E=(\S+) ([0-9.]+), ([^,]+), ([^,]+), (\S+)")


def points_du_journal(kid):
    """{(S, E) : {gamma : (Wp_c, k_c, Wi_c)}} lus dans les journaux du kernel."""
    out = {}
    for f in glob.glob(os.path.join(KO, kid, "*.log")):
        txt = open(f, encoding="utf-8", errors="replace").read()
        try:                                   # journal Kaggle = liste JSON
            flux = " ".join(e.get("data", "") for e in json.loads(txt))
        except Exception:
            flux = txt
        for S, E, g, w, k, wi in ROW.findall(flux):
            cle = (float(S), float(E))
            out.setdefault(cle, {})[round(float(g), 5)] = (w.strip(), k.strip(), wi.strip())
    return out


def deja_dans_les_csv(S, E):
    vus = set()
    for f in glob.glob(os.path.join(RG, "tail_S%.2f_E%g_p*.csv" % (S, E))):
        for line in open(f).read().splitlines()[1:]:
            if line.strip():
                vus.add(round(float(line.split(",")[0]), 5))
    return vus


def main():
    total = {}
    for d in sorted(glob.glob(os.path.join(KO, "shqgt-*"))):
        for cle, pts in points_du_journal(os.path.basename(d)).items():
            total.setdefault(cle, {}).update(pts)
    for (S, E), pts in sorted(total.items()):
        vus = deja_dans_les_csv(S, E)
        neufs = {g: v for g, v in pts.items() if g not in vus}
        print("S=%-5g E=%-5g : %3d points dans les journaux, %3d absents des CSV"
              % (S, E, len(pts), len(neufs)))
        if not neufs:
            continue
        f = os.path.join(RG, "tail_S%.2f_E%g_plog.csv" % (S, E))
        with open(f, "w") as fh:
            fh.write("gamma, Wp_c, alpha_c, Wi_c_eq\n")
            for g in sorted(neufs):
                w, k, wi = neufs[g]
                fh.write("%.5f, %s, %s, %s\n" % (g, w, k, wi))
        print("   ecrit", f)


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""Met dans l'article les figures des cas TERMINES seulement.

Demande de l'auteur (18/09) : ne pas afficher une figure tant que son cas n'est
pas calculé à 100 %. Un cas est complet quand resultats_G/G_S<S>_E<E>.csv porte
ses 19 points de basse fréquence, ses 375 points de corps et ses 62 points de
queue.

Pour chacun des quatre blocs figure du maître (figures 2 à 5, un E par bloc,
quatre panneaux S = 1, 0.9, 0.7, 0.3), ce script écrit :
  * \\includegraphics{figures/G_E<E>_S<S>.png} si le cas est complet,
  * le cadre vide « computation in progress » sinon.
Il est idempotent et se relance à chaque arrivée de données :
    python _G_figs_article.py && python _G_inserer_figures.py
puis, dans article_JFM/ : python _build_prf.py, pdflatex.

Usage : python _G_inserer_figures.py
"""
import os
import re

import numpy as np

RG = "resultats_G"
MAITRE = os.path.join("article_JFM", "article_JFM_full.tex")
SS = (0.98, 0.9, 0.7, 0.3)
BLOCS = (("fig:E0p01", 0.01), ("fig:E1", 1.0), ("fig:E10", 10.0), ("fig:E50", 50.0))
VIDE = (r"\fbox{\parbox[c][0.62\textwidth][c]{0.92\textwidth}"
        r"{\centering\small computation in progress}}")
PANNEAU = re.compile(r"\\includegraphics\[width=\\textwidth\]\{figures/G_E[^}]*\}"
                     r"|" + re.escape(VIDE))
N_BF, N_CORPS, N_QUEUE = 19, 375, 62


def tag(v):
    return ("%g" % v).replace(".", "p")


def complet(S, E):
    f = os.path.join(RG, "G_S%.2f_E%g.csv" % (S, E))
    if not os.path.exists(f):
        return False, (0, 0, 0)
    d = np.atleast_2d(np.genfromtxt(f, delimiter=",", skip_header=1))
    g = d[np.isfinite(d[:, 1])][:, 0]
    n = (int((g < 0.1).sum()), int(((g >= 0.1) & (g < 7.78)).sum()), int((g >= 7.78).sum()))
    return n >= (N_BF, N_CORPS, N_QUEUE) and n[0] >= N_BF and n[1] >= N_CORPS and n[2] >= N_QUEUE, n


def main():
    t = open(MAITRE, encoding="utf-8").read()
    change = 0
    for lab, E in BLOCS:
        i = t.index("\\label{%s}" % lab)
        j = t.rindex("\\begin{figure}", 0, i)
        bloc = t[j:i]
        trouves = PANNEAU.findall(bloc)
        assert len(trouves) == 4, (lab, len(trouves))
        neufs = []
        for S in SS:
            ok, n = complet(S, E)
            neufs.append(r"\includegraphics[width=\textwidth]{figures/G_E%s_S%s.png}"
                         % (tag(E), tag(S)) if ok else VIDE)
            etat = "figure" if ok else "cadre vide (bf %d/%d, corps %d/%d, queues %d/%d)" % (
                n[0], N_BF, n[1], N_CORPS, n[2], N_QUEUE)
            print("E=%-5g S=%-5g : %s" % (E, S, etat))
        it = iter(neufs)
        neuf_bloc = PANNEAU.sub(lambda m: next(it), bloc)
        if neuf_bloc != bloc:
            change += 1
        t = t[:j] + neuf_bloc + t[i:]
    open(MAITRE, "w", encoding="utf-8").write(t)
    print("\nblocs modifies :", change, "; cadres vides restants :",
          t.count("computation in progress"))


if __name__ == "__main__":
    main()

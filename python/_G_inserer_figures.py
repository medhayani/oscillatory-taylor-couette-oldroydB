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
# Cadre vide a la hauteur d'un panneau empile (disposition 4 x 1 du 20/09) ;
# l'ancien cadre, haut de 0.62 textwidth, datait de la disposition 2 x 2.
VIDE = (r"\fbox{\parbox[c][0.24\textwidth][c]{0.92\textwidth}"
        r"{\centering\small computation in progress}}")
PANNEAU = re.compile(r"\\includegraphics\[width=\\textwidth\]\{figures/G_E[^}]*\}"
                     r"|\\fbox\{\\parbox\[c\]\[0\.\d+\\textwidth\]\[c\]\{0\.92\\textwidth\}"
                     r"\{\\centering\\small computation in progress\}\}")
N_BF, N_CORPS, N_QUEUE = 19, 375, 62


def tag(v):
    return ("%g" % v).replace(".", "p")


def complet(S, E):
    """Un cas est complet si _G_figs_article.py a juge sa courbe complete sur
    la plage tracee : ce script-la ecrit la figure, ou la retire sinon. Le
    critere vit donc a un seul endroit (20/09 : donnees assemblees GH)."""
    f = os.path.join("article_JFM", "figures", "G_E%s_S%s.png" % (tag(E), tag(S)))
    return os.path.exists(f), (0, 0, 0)


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
            etat = "figure" if ok else "cadre vide (courbe incomplete sur la plage tracee)"
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

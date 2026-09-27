# -*- coding: utf-8 -*-
"""Tableau des constantes de la loi basse frequence (annexe C.3).

A(S) = lim_{gamma->0} Wp_c gamma^2 pour le systeme SANS l'inertie de la
perturbation, et le nombre d'onde k_c qui realise ce minimum. La ligne rouge
en tirets des figures 2 a 5 est Wi_c = A(S) / (sqrt(eps) gamma^2).

Source : resultats_F/lowg_S<S>_E<E>.csv (campagne F', 06/09/2026), au plus
petit gamma de chaque fichier, moyennee sur les quatre E.

ATTENTION (23/09/2026). Ce script lisait la RACINE du dossier de travail, ou
traine une serie lowg_*.csv du 21/08 qui donne A = 7.54 / 8.03 / 9.55 / 16.11,
soit un facteur 1.37 au-dessus. Ce n'est pas la serie de l'article : le relancer
tel quel aurait change les constantes du manuscrit sans prevenir. D'ou le
chemin explicite et l'assertion sur les valeurs attendues.

Sortie : article_JFM/wkb_table.tex
Usage : python _gen_wkb_table.py
"""
import os

import numpy as np

SRC = "resultats_F"
SS = (0.98, 0.9, 0.7, 0.3)
LAB = {0.98: "1", 0.9: "0.9", 0.7: "0.7", 0.3: "0.3"}
ES = (0.01, 1.0, 10.0, 50.0)
A_ATTENDU = {0.98: 5.50, 0.9: 5.88, 0.7: 6.99, 0.3: 11.47}   # garde-fou


def main():
    lignes = ["%% genere par _gen_wkb_table.py depuis %s (campagne F')" % SRC,
              r"\begin{tabular}{ccc}", r"\toprule",
              r"$S$ & $A$ & $k_c$\\", r"\midrule"]
    print("S     A        dispersion sur E (%)   k_c")
    for S in SS:
        A, K = [], []
        for E in ES:
            f = os.path.join(SRC, "lowg_S%.2f_E%g.csv" % (S, E))
            if not os.path.exists(f):
                continue
            d = np.atleast_2d(np.genfromtxt(f, delimiter=",", skip_header=1))
            A.append(float(d[0, 1] * d[0, 0] ** 2))
            K.append(float(d[0, 2]))
        assert A, "aucun fichier lowg dans %s pour S = %g" % (SRC, S)
        a = float(np.mean(A))
        disp = 100.0 * (max(A) - min(A)) / a
        assert abs(a - A_ATTENDU[S]) < 0.02, (
            "A = %.4f a S = %g : ce n'est pas la valeur du manuscrit (%.2f). "
            "Verifier que %s est bien la campagne F'." % (a, S, A_ATTENDU[S], SRC))
        assert len(set(K)) == 1, "k_c varie avec E a S = %g : %s" % (S, K)
        lignes.append(r"$%s$ & $%.2f$ & $%.1f$\\" % (LAB[S], a, K[0]))
        print("%-5s %-8.4f %-22.4f %.1f" % (LAB[S], a, disp, K[0]))
    lignes += [r"\bottomrule", r"\end{tabular}"]
    out = os.path.join("article_JFM", "wkb_table.tex")
    open(out, "w", encoding="utf-8").write("\n".join(lignes) + "\n")
    print("ecrit", out)


if __name__ == "__main__":
    main()

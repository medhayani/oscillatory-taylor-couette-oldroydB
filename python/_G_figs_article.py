# -*- coding: utf-8 -*-
"""Figures de la campagne G : UNE figure par cas (S, E), sans titre.

Rendu LaTeX reel (newtxtext/newtxmath, comme les figures de la campagne F'),
axes et etiquettes en grand. Wi_c(gamma) en noir (log-log, axe de gauche),
k_c(gamma) en bleu (axe de droite, lineaire). Pas de marque au minimum
(demande de l'auteur, 19/09) : les valeurs sont dans minima.txt.
k_c est trace en tirets la ou il touche le bord de la grille en k (40 sous
gamma = 7.8, 80 au-dessus ; 80 aussi pour les points rebalayes, E = 0.01 et
gamma <= 0.40) : le minimum en k n'y est pas atteint et Wi_c y est un majorant.
Les courbes a E = 0.01 commencent a gamma = 0.18 (G_MIN) : en dessous,
l'entrefer etroit sort de son domaine et N = 28 n'est pas converge. Pas de bande grise, pas de loi WKB en gamma^-2 (annexe C.3 : elle
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
from matplotlib.ticker import LogLocator, NullLocator, NullFormatter

RG = "resultats_G"
OUT = os.path.join(RG, "figures_G")
FIG_ART = os.path.join("article_JFM", "figures")
SS = (0.98, 0.9, 0.7, 0.3)
ES = (0.01, 1.0, 10.0, 50.0)
C_WI, C_K, GRID, C_WKB = "#000000", "#1f5fa9", "#e3e2dc", "#c0392b"
# Loi sans inertie (annexe C.3) : Wi_c = eps^{-1/2} A(S) / gamma^2, avec A(S)
# de article_JFM/wkb_table.tex (campagne quasi statique). Tracee en tirets
# rouges de gamma = 0.1 a 0.8 (demande de l'auteur, 19/09), comme repere : le
# systeme complet a E fixe ne la suit pas (pente mesuree -1, non -2).
A_WKB = {0.98: 5.50, 0.9: 5.88, 0.7: 6.99, 0.3: 11.47}
X0, X1 = 0.1, 22.0          # figures a partir de gamma = 0.1 (demande de l'auteur, 19/09)
GAP = 1.35
LOG = True           # True : axes logarithmiques ; False : axes lineaires
KMAX_BAS, KMAX_QUEUE = 40.0, 80.0
SAUT_MAX = 2.5       # facteur maximal tolere entre deux points voisins d'une courbe

# Debut de la COURBE, par elasticite (l'axe, lui, part toujours de X0 = 0.1).
# A E = 0.01 les points sous gamma = 0.18 ne sont pas utilisables : d'une part
# gamma^2 y devient comparable a eps^2 (l'entrefer etroit sort de son domaine),
# d'autre part la collocation a N = 28 du balayage n'y est pas convergee --
# verifie le 23/09 a S = 0.7, k = 57 : a gamma = 0.14, sigma reste negatif de
# Wp = 400 a 3600 avec N = 28, alors qu'il est POSITIF a Wp = 400 avec N = 44 ;
# a gamma = 0.16, en revanche, le seuil vaut 102.3 a N = 28, 40, 44, 56 et 60.
G_MIN = {0.01: 0.18}

# Prolongement sous le premier point calcule, jusqu'a la valeur donnee ici.
# Ce n'est PAS un calcul : c'est la loi de puissance locale, de pente prise sur
# les NFIT premiers points reels et ANCREE sur le premier d'entre eux, donc
# continue en valeur et en pente -- le raccord ne se voit pas, comme demande
# par l'auteur le 23/09 (<< tu mets pas pointille, tu mets qu'il est lie avec
# la courbe >>). Wi_c et k_c sont prolonges tous les deux.
#
# Pourquoi une loi de puissance et non une courbure. Essaye et rejete : un
# ajustement quadratique en log-log donne, a gamma = 0.1, 506 / 397 / 658 /
# 1433 selon S -- il CASSE l'ordre Wi(0.3) > Wi(0.7) > Wi(0.9) > Wi(1), et il
# s'ecarte davantage des seuls points vraiment converges. La loi de puissance
# respecte l'ordre et passe a +0.1 % et -1.6 % des deux points converges a
# N = 56 (S = 0.7, k = 57 : Wi_c = 273 a gamma = 0.16 et 360 a gamma = 0.14).
# La courbe EST une droite en log-log dans cette plage ; ce n'est pas un defaut
# du trace. Le prolongement reste annonce dans la legende et dans le texte.
EXTRAP = {0.01: 0.1}
NFIT = 3
KLIM_EXT = (0.0, 132.0)          # axe k_c quand k_c est prolonge (il monte a 120 a S = 0.3)


def prolonge(x, y, cible, nfit=NFIT, log_y=True):
    """Prolonge (x, y) jusqu'a cible par la pente locale des nfit premiers
    points, ancree sur le premier : continue en valeur et en pente."""
    lx = np.log(x[:nfit])
    ly = np.log(y[:nfit]) if log_y else np.asarray(y[:nfit], float)
    pente = np.polyfit(lx, ly, 1)[0]
    xe = np.geomspace(cible, x[0], 30)
    ye = (ly[0] if not log_y else np.log(y[0])) + pente * (np.log(xe) - np.log(x[0]))
    return xe, (np.exp(ye) if log_y else ye)


# Lissage d'AFFICHAGE, sur [10, 20] seulement (demande de l'auteur, 27/09 :
# << un filtre a partir de gamma = 10 jusqu'a 20, tu mets pas les parasites >>).
# Mediane glissante sur FEN_MED points -- elle ote les pics isoles et garde les
# marches -- puis moyenne glissante sur FEN_MOY points. Les fenetres debordent
# sur les voisins bruts, donc le raccord a gamma = 10 est continu. Les DONNEES
# ne sont pas touchees : minima.txt, le tableau 3 et le texte lisent les
# fichiers bruts. Declare dans la legende de la figure 2.
LISSAGE = (10.0, 20.0)
# 5 + 3 laissait 8 a 10 sauts de k_c > 10 sur les queues S = 0.98 a E >= 1
# (alternance de branches sur 2 a 4 points de suite) ; 11 + 5 n'en laisse
# aucun, ni aucun saut de log Wi_c > 0.15, sur les 16 cas (mesure du 27/09),
# avec une fenetre qui reste locale : +-1.0 en gamma pour la mediane.
FEN_MED, FEN_MOY = 11, 5


def lisse(g, y, log_y=False):
    """Copie de y lissee la ou g est dans LISSAGE ; ailleurs y est rendu tel
    quel. Le lissage de Wi_c se fait en logarithme (log_y=True)."""
    y = np.asarray(y, float)
    o = np.argsort(g)
    gs = g[o]
    v = np.log(y[o]) if log_y else y[o].copy()
    n = len(v)
    lo, hi = LISSAGE
    idx = np.flatnonzero((gs >= lo - 1e-9) & (gs <= hi + 1e-9))
    if idx.size < FEN_MED:
        return y
    med = v.copy()
    h = FEN_MED // 2
    for i in idx:
        med[i] = np.nanmedian(v[max(0, i - h):min(n, i + h + 1)])
    out = med.copy()
    h2 = FEN_MOY // 2
    for i in idx:
        out[i] = np.nanmean(med[max(0, i - h2):min(n, i + h2 + 1)])
    res = np.exp(out) if log_y else out
    back = np.empty_like(res)
    back[o] = res
    return back


STYLE_LATEX = {
    "text.usetex": True,
    "font.family": "serif",
    "text.latex.preamble": r"\usepackage{newtxtext,newtxmath}",
    "font.size": 16, "axes.labelsize": 22,
    "xtick.labelsize": 17, "ytick.labelsize": 17,
    "axes.linewidth": 1.0, "axes.edgecolor": "#222222",
    "xtick.major.width": 1.0, "ytick.major.width": 1.0,
    "xtick.major.size": 5, "ytick.major.size": 5,
}
STYLE_SANS_LATEX = dict(STYLE_LATEX, **{"text.usetex": False,
                                        "mathtext.fontset": "stix"})


def tag(v):
    return ("%g" % v).replace(".", "p")


def load(S, E):
    """Courbe assemblee GH (campagne H corrigee sous De = 0.3, G au-dessus,
    ecrite par _GH_merge.py) ; a defaut, la campagne G seule."""
    f = os.path.join(RG, "GH_S%.2f_E%g.csv" % (S, E))
    if not os.path.exists(f):
        f = os.path.join(RG, "G_S%.2f_E%g.csv" % (S, E))
    if not os.path.exists(f):
        return None
    d = np.atleast_2d(np.genfromtxt(f, delimiter=",", skip_header=1))
    d = d[np.isfinite(d[:, 1])]
    d = d[d[:, 0] >= G_MIN.get(E, X0) - 1e-9]
    return (d[:, 0], d[:, 3], d[:, 2]) if len(d) else None


def complet(d, E):
    """Cas complet SUR LA PLAGE TRACEE, de G_MIN(E) a 20 : il commence la,
    finit a 20, et n'a aucun trou (pas de 0.02 dans le corps, 0.2 en queue ;
    les grilles ont des sauts de 0.04 aux limites de sous-blocs)."""
    g0 = G_MIN.get(E, X0)
    g = np.sort(d[0][d[0] >= g0 - 1e-9])
    if len(g) < 50 or g[0] > g0 + 1e-3 or g[-1] < 19.99:
        return False
    pas = np.diff(g)
    corps = pas[g[1:] < 7.79]
    queue = pas[g[1:] >= 7.79]
    if not (corps.max() <= 0.0601 and (len(queue) == 0 or queue.max() <= 0.2501)):
        return False
    # Presentable : pas de saut d'un facteur > SAUT_MAX entre deux points voisins.
    # Un tel saut n'est pas un seuil, c'est la fenetre en k qui a change de
    # branche (regle de Codex : ne pas lisser) ; le balayage complet en k
    # (_gen_kscan_kaggle.py) remplace ces points.
    o = np.argsort(d[0])
    w = d[1][o][d[0][o] >= g0 - 1e-9]
    r = w[1:] / w[:-1]
    return bool(np.max(np.abs(np.log(r))) <= np.log(SAUT_MAX))


def kmax_of(g, E=None):
    """Bord de la grille en k, selon la famille dont vient le point.

    La campagne balaie k <= 40 sous gamma = 7.78 et k <= 80 au-dessus. Mais a
    E = 0.01 les points gamma <= 0.40 ont ete REBALAYES en entier, k = 1 a 80
    pas 1 (kernels shqk2-*, 23/09), avec un minimum interieur verifie : leur
    bord est donc 80, pas 40, et k_c y est un vrai minimum."""
    km = np.where(g < 7.78, KMAX_BAS, KMAX_QUEUE)
    if E is not None and E < 0.1:
        km = np.where(g <= 0.405, KMAX_QUEUE, km)
    return km


# Traits verticaux aux changements de mode : RETIRES PAR L'AUTEUR le 23/09
# (<< tu enleves les lignes >>). La detection est conservee ci-dessous ; mettre
# MODES = True pour les retracer. Ne pas les remettre sans demande explicite.
MODES = False
C_MODE = "#8a8a8a"       # trait vertical des changements de mode
M_FEN = 5                # demi-fenetre (en points) de la mediane et des pentes
M_DK = 4.0               # saut de k_c au-dela duquel on parle de changement
M_COUDE = 0.35           # coude minimal de la pente de log Wi_c en log gamma


def transitions(g, wi, k, E):
    """Frequences ou la courbe passe d'un mode a un autre.

    Un changement de mode se lit deux fois : le nombre d'onde selectionne saute
    d'une branche a l'autre, et Wi_c fait un petit coude. On exige les deux, et
    on exige que le saut de k_c soit DURABLE -- mediane sur M_FEN points avant
    et apres -- sinon le va-et-vient de k_c au bord de la grille (cas S = 1 a
    E >= 1) donnerait des dizaines de fausses marques. Les zones ou k_c est au
    bord de la grille sont ecartees : le minimum en k n'y est pas atteint.
    """
    sat = k >= kmax_of(g, E) - 1e-9
    lw, lg = np.log(wi), np.log(g)
    m = M_FEN
    flags = []
    for i in range(m, len(g) - m - 1):
        if sat[i - m:i + m + 2].any():
            continue
        if abs(np.median(k[i + 1:i + 1 + m]) - np.median(k[i - m:i])) < M_DK:
            continue
        pa = (lw[i] - lw[i - m]) / (lg[i] - lg[i - m])
        pb = (lw[i + 1 + m] - lw[i + 1]) / (lg[i + 1 + m] - lg[i + 1])
        if abs(pb - pa) < M_COUDE:
            continue
        flags.append(i)
    grp = []                              # fusionner les indices voisins
    for i in flags:
        if grp and i - grp[-1][-1] <= m:
            grp[-1].append(i)
        else:
            grp.append([i])
    return [float(np.exp(np.mean(lg[gg]))) for gg in grp]


def cut(g, *ys):
    """Coupe le trait la ou la couverture en gamma a un trou."""
    if len(g) < 2:
        return (g,) + ys
    c = np.flatnonzero(g[1:] / g[:-1] > GAP) + 1
    return (np.insert(g, c, np.nan),) + tuple(np.insert(y, c, np.nan) for y in ys)


# Geometrie des panneaux empiles, en pouces. Le cadre a exactement la meme
# taille sur les quatre panneaux ; seul le dernier porte les graduations et le
# titre de l'axe gamma, d'ou une image un peu plus haute pour lui. Les marges
# sont donc reduites au minimum entre deux panneaux.
L_FIG, L_GAUCHE, L_CADRE = 7.6, 0.80, 6.00      # largeurs
H_CADRE, H_TITRE, H_BAS, H_BAS_DERNIER = 1.70, 0.30, 0.08, 0.62


def dessine(d, chemin, S=None, E=None):
    g, wi, k = d
    dernier = (S == 0.3)
    h_bas = H_BAS_DERNIER if dernier else H_BAS
    h_fig = H_CADRE + H_TITRE + h_bas
    fig = plt.figure(figsize=(L_FIG, h_fig), facecolor="white")
    ax = fig.add_axes([L_GAUCHE / L_FIG, h_bas / h_fig,
                       L_CADRE / L_FIG, H_CADRE / h_fig])
    if LOG:
        ax.set_xscale("log")
        ax.set_yscale("log")
        # Le cadre commence a la PREMIERE valeur de gamma du cas, comme pour
        # les autres elasticites : pas de bande vide a gauche (demande de
        # l'auteur, 23/09). Les quatre panneaux d'une meme figure partagent
        # cette borne, puisqu'ils partagent G_MIN.
        ax.set_xlim(EXTRAP.get(E, G_MIN.get(E, X0)), X1)
        # Graduations de l'axe gamma : 0.1, 0.5, 1, 5, 10, 20 en clair (demande
        # de l'auteur, 27/09) ; les graduations mineures gardent leurs traits
        # mais pas d'etiquette.
        ax.set_xticks([0.1, 0.5, 1.0, 5.0, 10.0, 20.0])
        ax.set_xticklabels(["0.1", "0.5", "1", "5", "10", "20"])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_ylim(float(np.nanmin(wi)) / 1.6, float(np.nanmax(wi)) * 1.6)   # cadre sur les donnees
    else:                                   # axes lineaires (demande du 17/09)
        ax.set_xlim(0.0, 20.5)
        ax.set_ylim(0.0, float(np.nanmax(wi)) * 1.05)
    ax.grid(True, which="major", color=GRID, lw=0.7)
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=4))     # peu de graduations : panneaux etroits
    ax.yaxis.set_minor_locator(NullLocator())
    if dernier:                        # graduations et titre gamma : dernier panneau seulement
        ax.set_xlabel(r"$\gamma$")
    else:
        ax.tick_params(axis="x", labelbottom=False)
    ETIQ = {0.98: ("a", "1"), 0.9: ("b", "0.9"), 0.7: ("c", "0.7"), 0.3: ("d", "0.3")}
    if S in ETIQ:                      # (a) S = ... dans le panneau, plus sous la figure
        lettre, val = ETIQ[S]
        ax.set_title(r"\textbf{(\textit{%s})\quad $S=%s$}" % (lettre, val),
                     loc="center", fontsize=17, pad=2)    # centre, en gras, au-dessus du cadre
    ax.set_ylabel(r"$Wi_c$")
    axk = ax.twinx()
    axk.set_ylabel(r"$k_c$", color=C_K)
    axk.tick_params(axis="y", colors=C_K)
    if E in EXTRAP:                    # k_c prolonge : il passe 80
        axk.set_ylim(*KLIM_EXT)
        axk.set_yticks([0, 60, 120])
    else:
        axk.set_ylim(0, 85)
        axk.set_yticks([0, 40, 80])

    if S in A_WKB and E is not None:
        # Loi sans inertie, tracee du bord gauche du cadre a gamma = 0.8 sur
        # CHAQUE cas (demande de l'auteur, 21/09, reprise le 23/09 : elle doit
        # commencer a 0.1 comme la courbe). Le cadre est ouvert vers le haut
        # pour qu'elle y entre en entier.
        g_lo, g_hi = max(EXTRAP.get(E, G_MIN.get(E, X0)), 0.1), 0.8
        if g_hi > g_lo:
            gw = np.geomspace(g_lo, g_hi, 40)
            ww_wkb = A_WKB[S] / np.sqrt(0.14) / gw ** 2
            ax.plot(gw, ww_wkb, ls=(0, (5, 3)),
                    lw=2.4, color=C_WKB, zorder=2)   # identifiee dans la legende de la figure
            if LOG:
                bas, haut = ax.get_ylim()
                ax.set_ylim(bas, max(haut, float(ww_wkb.max()) * 1.25))
    if MODES:                                          # changements de mode
        o0 = np.argsort(g)
        for gt in transitions(g[o0], wi[o0], k[o0], E):
            ax.axvline(gt, color=C_MODE, lw=0.9, zorder=1)
    wi_l, k_l = lisse(g, wi, log_y=True), lisse(g, k)     # lissage d'affichage sur [10, 20]
    gg, ww, kk = cut(g, wi_l, k_l)
    ax.plot(gg, ww, "-", color=C_WI, lw=2.0, zorder=4)
    axk.plot(gg, kk, "-", color=C_K, lw=1.2, zorder=3)
    if E in EXTRAP and len(g) > NFIT:
        # Prolongement, dans le meme trait que la courbe (il s'y raccorde en
        # valeur et en pente). Annonce dans la legende de la figure.
        o = np.argsort(g)
        gs, ws, ks = g[o], wi[o], k[o]
        ge, we = prolonge(gs, ws, EXTRAP[E])
        _, ke = prolonge(gs, ks, EXTRAP[E], log_y=False)
        ax.plot(ge, we, "-", color=C_WI, lw=2.0, zorder=4)
        axk.plot(ge, ke, "-", color=C_K, lw=1.2, zorder=3)
        if LOG:                      # le cadre doit contenir le prolongement
            bas, haut = ax.get_ylim()
            ax.set_ylim(bas, max(haut, float(we.max()) * 1.6))
    sat = k >= kmax_of(g, E) - 1e-9
    if sat.any():
        gs, kss = cut(g, np.where(sat, k_l, np.nan))   # trait sur la courbe lissee, detection sur k brut
        axk.plot(gs, kss, ls=(0, (4, 2)), color=C_K, lw=2.0, zorder=5)
    i = int(np.argmin(wi))          # garde pour le tableau des minima, plus marque sur la figure
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
            if not complet(d, E):     # cas non calcule a 100 % : pas de figure
                for f in (os.path.join(OUT, nom), os.path.join(FIG_ART, nom)):
                    if os.path.exists(f):
                        os.remove(f)
                g = np.sort(d[0][d[0] >= G_MIN.get(E, X0) - 1e-9])
                print("E=%-5g S=%-5g : incomplet sur [%g, 20] (premier gamma trace %.3g, "
                      "plus grand trou %.3g), figure retiree"
                      % (E, S, G_MIN.get(E, X0), g[0] if len(g) else np.nan,
                         np.diff(g).max() if len(g) > 1 else np.nan))
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

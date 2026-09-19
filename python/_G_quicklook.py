# -*- coding: utf-8 -*-
"""Campagne G : fusion des blocs recus et comparaison a la campagne F'.

Sources G (resultats_G/) : lowg_*_p*.csv (gamma 0.010-0.089), fine_*_p*.csv
(0.10-7.75), tail_*_p*.csv (7.8-20) -> consolide resultats_G/G_S<S>_E<E>.csv
(gamma, Wp_c, k_c, Wi_c).
Reference F' : resultats_F/fine_S<S>_E<E>.csv (0.01-20, basse frequence et
queues incluses).
Figures : resultats_G/G_vs_F_S<S>.png, une par S ; colonnes E = 0.01, 1, 10, 50 ;
lignes Wi_c(gamma), k_c(gamma), rapport Wi_c(G) / Wi_c(F').
Resume des minima : resultats_G/G_minima.txt.
Palette validee (dataviz/validate_palette.js, mode clair : ALL PASS) :
G bleu #2a78d6, F' orange #eb6834 ; encres et grille neutres.
Relancable a chaque pull. Usage : python _G_quicklook.py
"""
import os, glob
import numpy as np

RG, RF = "resultats_G", "resultats_F"
SQ = np.sqrt(0.14)
SS = (0.98, 0.9, 0.7, 0.3)
ES = (0.01, 1.0, 10.0, 50.0)
C_G, C_F = "#2a78d6", "#eb6834"
SURF, INK1, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
GAP = 1.35          # rupture de trait si deux gamma consecutifs sont dans un rapport > GAP


def read_parts(S, E):
    rows, fams = {}, {}
    for fam in ("lowg", "fine", "tail"):
        n = 0
        for f in sorted(glob.glob(os.path.join(RG, "%s_S%.2f_E%g_p*.csv" % (fam, S, E)))):
            for line in open(f).read().splitlines()[1:]:
                v = [x.strip() for x in line.split(",")]
                if len(v) < 3 or not v[0]:
                    continue
                g = round(float(v[0]), 5)
                w = float(v[1]) if v[1] != "nan" else np.nan
                k = float(v[2]) if v[2] != "nan" else np.nan
                rows[g] = (w, k)
                n += 1
        fams[fam] = n
    if not rows:
        return None, fams
    g = np.array(sorted(rows))
    w = np.array([rows[x][0] for x in g])
    k = np.array([rows[x][1] for x in g])
    with open(os.path.join(RG, "G_S%.2f_E%g.csv" % (S, E)), "w") as fh:
        fh.write("gamma,Wp_c,k_c,Wi_c\n")
        for a, b, c in zip(g, w, k):
            fh.write("%.5f,%.6g,%g,%.6g\n" % (a, b, c, b / SQ))
    return (g, w, k), fams


def load_F(S, E):
    f = os.path.join(RF, "fine_S%.2f_E%g.csv" % (S, E))
    if not os.path.exists(f):
        return None
    d = np.atleast_2d(np.genfromtxt(f, delimiter=",", skip_header=1))
    return d[:, 0], d[:, 1], d[:, 2]


def with_gaps(g, *ys):
    """Insere des NaN aux trous de couverture pour ne pas relier deux familles."""
    if len(g) < 2:
        return (g,) + ys
    cut = np.flatnonzero(g[1:] / g[:-1] > GAP) + 1
    g2 = np.insert(g.astype(float), cut, np.nan)
    return (g2,) + tuple(np.insert(y.astype(float), cut, np.nan) for y in ys)


def minimum(g, w, k):
    ok = np.isfinite(w)
    if not ok.any():
        return None
    i = np.flatnonzero(ok)[int(np.argmin(w[ok]))]
    return g[i], w[i] / SQ, k[i]


def ratio(G, F):
    g, w, _ = G
    gf, wf, _ = F
    okf = np.isfinite(wf) & (wf > 0)
    if okf.sum() < 2:
        return None
    m = np.isfinite(w) & (g >= gf[okf].min()) & (g <= gf[okf].max())
    if not m.any():
        return None
    wi = np.exp(np.interp(np.log(g[m]), np.log(gf[okf]), np.log(wf[okf])))
    return g[m], w[m] / wi


def cover(fams):
    lab = (("lowg", "0.01-0.09"), ("fine", "0.1-7.75"), ("tail", "7.8-20"))
    return "   ".join("%s : %s" % (r, ("%d pts" % fams[f]) if fams.get(f) else "en cours") for f, r in lab)


def figure_S(S, data):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({
        "font.size": 8.5, "font.family": "DejaVu Sans",
        "axes.edgecolor": AXIS, "axes.linewidth": 0.7, "axes.labelcolor": INK2,
        "xtick.color": INK2, "ytick.color": INK2, "text.color": INK1,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
        "lines.solid_capstyle": "round", "lines.solid_joinstyle": "round"})
    fig, axes = plt.subplots(3, 4, figsize=(14.0, 8.6), sharex=True, facecolor=SURF,
                             gridspec_kw={"height_ratios": [1.35, 1.0, 0.85], "hspace": 0.12, "wspace": 0.26})
    for j, E in enumerate(ES):
        a1, a2, a3 = axes[:, j]
        for ax in (a1, a2, a3):
            ax.set_facecolor(SURF)
            ax.set_xscale("log")
            ax.set_xlim(0.008, 25)
            ax.grid(True, which="major", color=GRID, lw=0.6)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
        a1.set_yscale("log")
        a3.set_yscale("log")
        a1.set_title("E = %g" % E, color=INK1, fontsize=10.5, pad=16)
        F = load_F(S, E)
        if F is not None:
            gF, wF, kF = with_gaps(F[0], F[1], F[2])
            a1.plot(gF, wF / SQ, color=C_F, lw=1.2, zorder=2)
            a2.plot(gF, kF, color=C_F, lw=1.2, zorder=2)
        a3.axhline(1.0, color=AXIS, lw=0.9, zorder=1)
        item = data.get(E)
        mF = minimum(*F) if F is not None else None
        if item is None:
            a3.text(0.5, 0.5, "campagne G : calcul en cours", transform=a3.transAxes,
                    ha="center", va="center", color=MUTED, fontsize=9)
            a1.text(0.0, 1.02, "G : en attente", transform=a1.transAxes, ha="left", va="bottom", color=MUTED, fontsize=7)
            if mF is not None:
                a1.plot([mF[0]], [mF[1]], "o", ms=5.5, mfc=C_F, mec=SURF, mew=1.2, zorder=4)
                a1.text(0.5, 0.97, "min F' : %.3g  (γ* = %.2f, k = %g)" % (mF[1], mF[0], mF[2]),
                        transform=a1.transAxes, ha="center", va="top", color=INK2, fontsize=7.4)
            continue
        G, fams = item
        g, w, k = G
        a1.text(0.0, 1.02, "G  " + cover(fams), transform=a1.transAxes, ha="left", va="bottom", color=MUTED, fontsize=6.6)
        gG, wG, kG = with_gaps(g, w, k)
        a1.plot(gG, wG / SQ, color=C_G, lw=1.5, zorder=3)
        a2.plot(gG, kG, color=C_G, lw=1.5, zorder=3)
        mG = minimum(g, w, k)
        box = []
        if mG is not None:
            a1.plot([mG[0]], [mG[1]], "o", ms=5.5, mfc=C_G, mec=SURF, mew=1.2, zorder=5)
            box.append("min G : %.3g  (γ* = %.2f, k = %g)" % (mG[1], mG[0], mG[2]))
        if mF is not None:
            a1.plot([mF[0]], [mF[1]], "o", ms=5.5, mfc=C_F, mec=SURF, mew=1.2, zorder=4)
            box.append("min F' : %.3g  (γ* = %.2f, k = %g)" % (mF[1], mF[0], mF[2]))
        if box:
            a1.text(0.5, 0.97, "\n".join(box), transform=a1.transAxes, ha="center", va="top",
                    color=INK1, fontsize=7.4, linespacing=1.5)
        if F is not None:
            r = ratio(G, F)
            if r is not None:
                gr, rr = with_gaps(r[0], r[1])
                a3.plot(gr, rr, color=C_G, lw=1.5, zorder=3)
                lo, hi = np.nanmin(rr), np.nanmax(rr)
                a3.set_ylim(min(0.5, lo / 1.3), max(2.0, hi * 1.3))
    axes[0, 0].set_ylabel("Wi_c", fontsize=9.5)
    axes[1, 0].set_ylabel("k_c", fontsize=9.5)
    axes[2, 0].set_ylabel("Wi_c(G) / Wi_c(F')", fontsize=9.5)
    for ax in axes[2, :]:
        ax.set_xlabel("γ", fontsize=10)
    name = "S = 0.98  (substitut de S = 1)" if S >= 0.98 else "S = %g" % S
    fig.suptitle(name, x=0.06, y=0.995, ha="left", fontsize=13, color=INK1, fontweight="bold")
    fig.legend(handles=[Line2D([0], [0], color=C_G, lw=2.0, label="inertie complète (campagne G)"),
                        Line2D([0], [0], color=C_F, lw=2.0, label="quasi statique (campagne F')")],
               loc="upper right", bbox_to_anchor=(0.985, 1.0), ncol=2, frameon=False, fontsize=9.5)
    fig.subplots_adjust(left=0.06, right=0.985, top=0.9, bottom=0.07)
    out = os.path.join(RG, "G_vs_F_S%.2f.png" % S)
    fig.savefig(out, dpi=130, facecolor=SURF)
    plt.close(fig)
    return out


def main():
    lines = ["S     E      points G (bf/corps/queue)  nan  | G : min Wi_c  gamma*   k_c* | F' : min Wi_c  gamma*   k_c* | G/F'"]
    figs = []
    for S in SS:
        data = {}
        for E in ES:
            G, fams = read_parts(S, E)
            if G is None:
                continue
            data[E] = (G, fams)
            g, w, k = G
            mG = minimum(g, w, k)
            F = load_F(S, E)
            mF = minimum(*F) if F is not None else None
            lines.append("%-5g %-5g  %4d/%4d/%4d             %4d  | %10.4g  %6.3f  %5.1f | %11.4g  %6.3f  %5.1f | %s"
                         % (S, E, fams.get("lowg", 0), fams.get("fine", 0), fams.get("tail", 0),
                            int((~np.isfinite(w)).sum()),
                            mG[1] if mG else np.nan, mG[0] if mG else np.nan, mG[2] if mG else np.nan,
                            mF[1] if mF else np.nan, mF[0] if mF else np.nan, mF[2] if mF else np.nan,
                            ("%.3f" % (mG[1] / mF[1])) if (mG and mF) else "-"))
        if data:
            figs.append(figure_S(S, data))
    txt = "\n".join(lines)
    open(os.path.join(RG, "G_minima.txt"), "w", encoding="utf-8").write(txt + "\n")
    print(txt)
    for f in figs:
        print("figure :", f)


if __name__ == "__main__":
    main()

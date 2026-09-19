# -*- coding: utf-8 -*-
"""Kernels Kaggle de la CAMPAGNE G : recalcul complet, gamma = 0.01 -> 20, avec
l'OPERATEUR COMPLET a huit champs (run_full_sweep.py) : inertie linearisee de
la perturbation retenue, coefficient C = Wp^2/(E S), termes -2 V v et u d_x V,
terme F' de (P6). Reference : article_JFM/version_PRFluids/
FORMULATION_EXACTE_apres_debat.tex.

Trois familles, sur les 16 cas S in {0.98, 0.9, 0.7, 0.3} x E in {0.01, 1, 10, 50} :
  * corps           shqg-s<S>e<E>-c1..3 : gamma 0.10 -> 7.75 pas 0.02, N = 28, k in [2, 40]
                    -> fine_S<S>_E<E>_p1..12.csv
  * basse frequence shqgl-s<S>e<E>      : 19 gamma log de 0.0100 a 0.0885, N = 28, k in [1, 40]
                    -> lowg_S<S>_E<E>_p1..4.csv        (ajout du 16/09)
  * queue           shqgt-s<S>e<E>-c1..2 : gamma 7.8 -> 20.0 pas 0.2, N = 48, k in [2, 80]
                    -> tail_S<S>_E<E>_p1..8.csv        (toutes les E, ajout du 16/09)
S = 1 exact (Sb = 0) donne des modes parasites (mesure du 16/09 : sigma brut
+363 a Wp = 0.05, queue de Chebyshev 0.12, aucun mode resolu) : S = 0.98, comme F'.
Amorcage : (k_c, Wp_c) de la campagne F' (resultats_F/) injectes (SEED) ; un
amorcage faux ne coute qu'un balayage log complet.
Ordre de poussee (ids.json) : corps, basse frequence, queues.

Isolation : dossier kg_G/, sorties rapatriees par _push_G.py dans resultats_G/.
Usage :  python _gen_G_kaggle.py
         python _push_G.py push | pull | list
"""
import os, json, shutil
import numpy as np

ROOT = "kg_G"
NSUB = 4                      # sous-blocs paralleles (4 vCPU par kernel)
SS = (0.98, 0.9, 0.7, 0.3)
ES = (0.01, 1.0, 10.0, 50.0)
CASES = [(S, E) for S in SS for E in ES]
# corps
G0, G1, DG = 0.10, 7.75, 0.02
NCH = 3
# basse frequence
LOW_G = np.round(np.geomspace(0.01, 0.1, 20)[:-1], 5)      # 19 points, 0.1 est dans le corps
# queues
TAIL_G = np.round(np.arange(7.8, 20.0 + 1e-9, 0.2), 2)     # 62 points
TNCH = 2
NCHEB_T, K_HI_T = 48, 80.0

FOOTER = '''

# ------------------------------------------------------------------ chunk
CASE_S, CASE_E = %s, %s
CH_G0, CH_G1 = %s, %s
NSUB = %d
OUTPAT = "%s"
SEED = %s


def _sub(t):
    g0, g1, k = t
    out = OUTPAT %% (CASE_S, CASE_E, k)
    dt, n = march(CASE_S, CASE_E, g0, g1, out)
    return "bloc %%d : %%d gammas en %%.0f s" %% (k, n, dt)


if __name__ == "__main__":
    import time
    from multiprocessing import Pool
    assert CORIOLIS and F_TERM, "operateur complet attendu (campagne G)"
    edges = np.linspace(CH_G0, CH_G1, NSUB + 1)
    tasks = []
    for k in range(NSUB):
        a = CH_G0 if k == 0 else round(edges[k] + DG, 3)
        tasks.append((round(a, 3), round(edges[k + 1], 3), %d + k))
    t0 = time.time()
    with Pool(NSUB) as pool:
        for msg in pool.imap_unordered(_sub, tasks):
            print(msg, flush=True)
    print("TERMINE chunk [%%.0f s]" %% (time.time() - t0))
    open("CHUNK_DONE.flag", "w").write("ok")
'''

FOOTER_LIST = '''

# ------------------------------------------------------------------ chunk (liste de gamma)
CASE_S, CASE_E = %s, %s
SUBLISTS = %s
BLOC0 = %d
OUTPAT = "%s"
SEED = %s


def _sub(t):
    gl, k = t
    out = OUTPAT %% (CASE_S, CASE_E, k)
    dt, n = march(CASE_S, CASE_E, None, None, out, gams=gl)
    return "bloc %%d : %%d gammas en %%.0f s" %% (k, n, dt)


if __name__ == "__main__":
    import time
    from multiprocessing import Pool
    assert CORIOLIS and F_TERM, "operateur complet attendu (campagne G)"
    tasks = [(gl, BLOC0 + i) for i, gl in enumerate(SUBLISTS)]
    t0 = time.time()
    with Pool(len(tasks)) as pool:
        for msg in pool.imap_unordered(_sub, tasks):
            print(msg, flush=True)
    print("TERMINE chunk [%%.0f s]" %% (time.time() - t0))
    open("CHUNK_DONE.flag", "w").write("ok")
'''


def tag(v):
    return ("%g" % v).replace(".", "p")


def build_body(ncheb=28, a_lo=2.0, a_hi=40.0):
    a = open("run_full_sweep.py", encoding="utf-8").read()
    assert "CORIOLIS = True" in a and "F_TERM = True" in a
    a = a[:a.index("def _selftest():")]
    for old, new in (("DG = 0.01\n", "DG = %g\n" % DG),
                     ("N = 28                  # points Chebyshev",
                      "N = %d                  # points Chebyshev" % ncheb),
                     ("A_LO, A_HI = 2.0, 40.0", "A_LO, A_HI = %g, %g" % (a_lo, a_hi))):
        assert old in a, old
        a = a.replace(old, new)
    return a


def seeds_for(S, E, starts, kind="fine"):
    """{gamma: (k_c, Wp_c)} de la campagne F' au plus pres de chaque depart."""
    f = os.path.join("resultats_F", "%s_S%.2f_E%g.csv" % (kind, S, E))
    out = {}
    if not os.path.exists(f):
        return out
    d = np.atleast_2d(np.genfromtxt(f, delimiter=",", skip_header=1))
    d = d[np.isfinite(d[:, 1])]
    if not len(d):
        return out
    for g in starts:
        i = int(np.argmin(np.abs(d[:, 0] - g)))
        out["%.5f" % g] = (float(d[i, 2]), float(d[i, 1]))
    return out


def write_kernel(kid, src):
    d = os.path.join(ROOT, kid)
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "script.py"), "w", encoding="utf-8", newline="\n").write(src)
    json.dump({
        "id": "USER/%s" % kid, "title": kid,
        "code_file": "script.py", "language": "python",
        "kernel_type": "script", "is_private": "true",
        "enable_gpu": "false", "enable_internet": "false",
        "dataset_sources": [], "competition_sources": [],
        "kernel_sources": [],
    }, open(os.path.join(d, "kernel-metadata.json"), "w"), indent=1)


def sub_starts(g0, g1):
    edges = np.linspace(g0, g1, NSUB + 1)
    return [round(g0 if k == 0 else edges[k] + DG, 3) for k in range(NSUB)]


def main():
    if os.path.isdir(ROOT):
        shutil.rmtree(ROOT)
    os.makedirs(ROOT)
    ids = []
    # ---- corps (noms et decoupage inchanges : les kernels deja pousses gardent leur etat)
    body = build_body()
    edges = np.linspace(G0, G1, NCH + 1)
    for S, E in CASES:
        for c in range(NCH):
            g0 = G0 if c == 0 else round(edges[c] + DG, 3)
            g1 = round(edges[c + 1], 3)
            kid = "shqg-s%se%s-c%d" % (tag(S), tag(E), c + 1)
            seed = {("%.5f" % float(k)): v for k, v in seeds_for(S, E, sub_starts(g0, g1)).items()}
            write_kernel(kid, body + FOOTER % (
                repr(S), repr(E), repr(g0), repr(g1), NSUB,
                "fine_S%.2f_E%g_p%d.csv", repr(seed), 1 + c * NSUB))
            ids.append(kid)
    # ---- basse frequence
    lbody = build_body(a_lo=1.0)
    for S, E in CASES:
        subs = [list(map(float, x)) for x in np.array_split(LOW_G, NSUB)]
        kid = "shqgl-s%se%s" % (tag(S), tag(E))
        seed = seeds_for(S, E, [s[0] for s in subs], kind="lowg")
        write_kernel(kid, lbody + FOOTER_LIST % (
            repr(S), repr(E), repr(subs), 1, "lowg_S%.2f_E%g_p%d.csv", repr(seed)))
        ids.append(kid)
    # ---- queues
    tbody = build_body(ncheb=NCHEB_T, a_hi=K_HI_T)
    for S, E in CASES:
        for c, part in enumerate(np.array_split(TAIL_G, TNCH)):
            subs = [list(map(float, x)) for x in np.array_split(part, NSUB)]
            kid = "shqgt-s%se%s-c%d" % (tag(S), tag(E), c + 1)
            seed = seeds_for(S, E, [s[0] for s in subs])
            write_kernel(kid, tbody + FOOTER_LIST % (
                repr(S), repr(E), repr(subs), 1 + c * NSUB, "tail_S%.2f_E%g_p%d.csv", repr(seed)))
            ids.append(kid)
    json.dump(ids, open(os.path.join(ROOT, "ids.json"), "w"), indent=1)
    nb = sum(k.startswith("shqg-") for k in ids)
    nl = sum(k.startswith("shqgl-") for k in ids)
    nt = sum(k.startswith("shqgt-") for k in ids)
    print("%d kernels dans %s/ : %d corps (0.10-7.75) + %d basse frequence (0.010-0.089) + %d queues (7.8-20)"
          % (len(ids), ROOT, nb, nl, nt))


if __name__ == "__main__":
    main()

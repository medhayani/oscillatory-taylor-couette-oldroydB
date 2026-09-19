# -*- coding: utf-8 -*-
"""Balayage du systeme COMPLET (P1)-(P8) a HUIT champs -- campagne G (16/09/2026).

Difference avec run_shaqfeh_sweep.py (campagne F') : l'inertie linearisee de
la perturbation est RETENUE dans (P1)-(P2), au coefficient exact
    C = eps^{1/2} Re Wp / S = Wp^2 / (E S)          (facteur We manquant dans l'article),
avec ses trois termes de meme ordre : d_t1, le couplage centrifuge -2 V v
(radial) et le transport du moment cinetique de base u d_x V (azimutal).
Reference : version_PRFluids/FORMULATION_EXACTE_apres_debat.tex, section 6.

    (P1) : C L2 d_t1 u = Sb Wp L2^2 u - k^2 Trr' - ik(dx^2+k^2) Trz
                         + k^2 Tzz' + k^2 Ttt  - 2 k^2 C V v
    (P2) : C d_t1 v    = Sb Wp L2 v + Trt' + ik Ttz  - C (d_x V) u
    (P3)-(P8) inchangees (terme F' de (P6) compris).

Etat y = [u, v, Trr, Trt, Trz, Ttt, Ttz, Tzz] (8N), conditions aux limites
u = u' = v = 0 aux parois eliminees exactement (y = Z q, 8N-6 inconnues) ;
M qdot = A(phi) q, propagation par exponentielle de matrice, m sous-pas.
S = 1 (Sb = 0) est REGULIER (plus de substitut S = 0.98).

Variables : Wp = eps^{1/2} We (controle), Sb = eta_s/eta_p = (1-S)/S,
E, gamma (De = 2 E gamma^2), alpha = k. Wi_c = eps^{-1/2} Wp_c.
Branche De < DE_SW : moyenne de cycle de la valeur propre dominante du
probleme generalise (A(phi), M) -- ecart mesure au Floquet exact 0.47*De
(_wkb_check.log), soit <= 2 % au raccord.

Sorties (memes conventions que la campagne F') :
  fine_S<S>_E<E>_p<k>.csv : gamma, Wp_c, alpha_c, Wi_c_eq
Usage :
  python run_full_sweep.py --test                     # auto-test
  python run_full_sweep.py <S> <E> <g0> <g1> <k>      # un bloc (Kaggle / local)
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
import sys, time, glob
import numpy as np
from numpy.linalg import solve, eig, lstsq
from scipy.linalg import expm
from scipy.optimize import brentq

CASE_S, CASE_E = None, None      # injectes par _gen_G_kaggle.py
SEED = {}                        # {"g0": (k_c, Wp_c)} injecte par le generateur
EPS = 0.14
N = 28                  # points Chebyshev
M_SUB = 40              # sous-pas expm / periode (m=40 vs 80 : <0.002 en sigma, verifie 15/09)
NSCAN = 22              # points du balayage log en Wp
WP_LO, WP_HI = 0.05, 3.0e3
DE_SW = 0.05            # De < DE_SW : branche adiabatique (moyenne de cycle)
NPH = 48                # phases de la moyenne de cycle
TAIL = 0.05             # filtre parasite : energie Chebyshev dans le tiers haut
CORIOLIS = True         # termes -2 V v et u d_x V (objet de la campagne G)
F_TERM = True           # terme -Wp^3 F' u de (P6)
EPSV = False            # terme eps V u de (M_theta), O(eps) : ecarte (entrefer etroit)


def chebdif(Npts):
    """D1, D2 sur x in [0,1] (Chebyshev-Gauss-Lobatto, x[0]=0, x[-1]=1)."""
    n = Npts - 1
    th = np.pi * np.arange(Npts) / n
    xi = np.cos(th)
    c = np.ones(Npts); c[0] = c[-1] = 2.0
    c *= (-1.0) ** np.arange(Npts)
    X = np.tile(xi, (Npts, 1)).T
    dX = X - X.T + np.eye(Npts)
    D = np.outer(c, 1.0 / c) / dX
    D -= np.diag(D.sum(axis=1))
    x = (1.0 - xi) / 2.0
    D1 = -2.0 * D
    return x, D1, D1 @ D1


def _tail_fraction(y, n):
    """Part d'energie de u dans le tiers superieur du spectre de Chebyshev."""
    xi = np.cos(np.pi * np.arange(n) / (n - 1))
    try:
        c = np.polynomial.chebyshev.chebfit(xi, y[:n][::-1], n - 1)
    except Exception:
        return 1.0
    e = np.abs(c) ** 2
    tot = e.sum()
    return float(e[2 * n // 3:].sum() / tot) if tot > 0 else 1.0


def _sigma_filtered(B, Z, n, tail=TAIL):
    """max Re(valeur propre) sur les modes bien resolus seulement."""
    try:
        w, V = eig(B)
    except np.linalg.LinAlgError:
        return 1e3
    for i in np.argsort(w.real)[::-1]:
        if not np.isfinite(w[i].real):
            continue
        if _tail_fraction(Z @ V[:, i], n) < tail:
            return float(w[i].real)
    return -1e3


class ShaqfehFull(object):
    """(P1)-(P8) a huit champs, inertie linearisee complete ; .sigma(Wp)."""

    def __init__(self, alpha, Sb, E, gam, n=None, m=None, steady=False):
        self.a, self.Sb, self.E, self.gam = float(alpha), float(Sb), float(E), float(gam)
        self.De = 2.0 * self.E * self.gam ** 2
        self.n = int(N if n is None else n)
        self.m = int(M_SUB if m is None else m)
        self.steady = bool(steady)
        n = self.n
        x, D1, D2 = chebdif(n)
        self.x, self.D1, self.D2 = x, D1, D2
        I = np.eye(n)
        self.L2 = D2 - self.a ** 2 * I
        Sb, De = self.Sb, self.De
        if steady:                       # Couette V = 1-x (LSM / Taylor)
            V1, V2 = 1.0 - x, np.zeros(n)
            V1p, V2p = -np.ones(n), np.zeros(n)
            S1, S2 = -np.ones(n), np.zeros(n)
            S1p, S2p = np.zeros(n), np.zeros(n)
            F0p, F1p, F2p = np.zeros(n), np.zeros(n), np.zeros(n)
        else:
            q2 = 2j * (1 + Sb) * gam ** 2 * (1 + 1j * De) / (1.0 + Sb * (1 + 1j * De))
            q = np.sqrt(q2)
            ch = np.cosh(q / 2.0)
            Vc = np.cosh(q * (x - 0.5)) / ch
            Vcp = q * np.sinh(q * (x - 0.5)) / ch
            Vcpp = q2 * Vc
            Sc = Vcp / (1 + 1j * De)
            Scp = Vcpp / (1 + 1j * De)
            V1, V2 = Vc.real, -Vc.imag
            V1p, V2p = Vcp.real, -Vcp.imag
            S1, S2 = Sc.real, -Sc.imag
            S1p, S2p = Scp.real, -Scp.imag
            if F_TERM:
                F0p = 2.0 * (Vcpp * np.conj(Vcp)).real / (1 + De ** 2)
                Fcp = 2.0 * Vcp * Vcpp / ((1 + 1j * De) * (1 + 2j * De))
                F1p, F2p = Fcp.real, -Fcp.imag
            else:
                F0p, F1p, F2p = np.zeros(n), np.zeros(n), np.zeros(n)
        self.V1, self.V2, self.V1p, self.V2p = V1, V2, V1p, V2p
        self.S1, self.S2, self.S1p, self.S2p = S1, S2, S1p, S2p
        self.F0p, self.F1p, self.F2p = F0p, F1p, F2p
        self._build_constraints()

    def C(self, Wp):
        """Coefficient inertiel reduit eps^{1/2} Re Wp / S = (1+Sb) Wp^2 / E."""
        return (1.0 + self.Sb) * Wp ** 2 / self.E

    def _build_constraints(self):
        n, N8 = self.n, 8 * self.n
        D1 = self.D1
        rows = np.zeros((6, N8))
        rows[0, 0] = 1.0; rows[1, n - 1] = 1.0
        rows[2, :n] = D1[0, :]; rows[3, :n] = D1[n - 1, :]
        rows[4, n] = 1.0; rows[5, 2 * n - 1] = 1.0
        dep = [0, 1, n - 2, n - 1, n, 2 * n - 1]
        ind = [i for i in range(N8) if i not in dep]
        W = -solve(rows[:, dep], rows[:, ind])
        Z = np.zeros((N8, len(ind)), dtype=complex)
        Z[ind, :] = np.eye(len(ind))
        Z[dep, :] = W
        self.Z = Z
        self.keep = ind

    def _blocks(self, Wp):
        """M, A0, Ac, As, A2c, A2s : M ydot = (A0 + Ac cos + As sin + A2c cos2 + A2s sin2) y."""
        n, a = self.n, self.a
        I = np.eye(n)
        D1, D2, L2 = self.D1, self.D2, self.L2
        Sb = self.Sb
        C = self.C(Wp)
        u, v = slice(0, n), slice(n, 2 * n)
        rr, rt, rz = slice(2 * n, 3 * n), slice(3 * n, 4 * n), slice(4 * n, 5 * n)
        tt, tz, zz = slice(5 * n, 6 * n), slice(6 * n, 7 * n), slice(7 * n, 8 * n)
        # ---- masse ----
        M = np.zeros((8 * n, 8 * n), dtype=complex)
        M[u, u] = C * L2
        M[v, v] = C * I
        for s in (rr, rt, rz, tt, tz, zz):
            M[s, s] = Wp * I
        # ---- partie constante ----
        A0 = np.zeros((8 * n, 8 * n), dtype=complex)
        A0[u, u] = Sb * Wp * (L2 @ L2)                    # (P1)
        A0[u, rr] = -a ** 2 * D1
        A0[u, rz] = -1j * a * (D2 + a ** 2 * I)
        A0[u, zz] = a ** 2 * D1
        A0[u, tt] = a ** 2 * I
        A0[v, v] = Sb * Wp * L2                           # (P2)
        A0[v, rt] = D1
        A0[v, tz] = 1j * a * I
        A0[rr, rr] = -I; A0[rr, u] = 2.0 * Wp * D1        # (P3)
        A0[rz, rz] = -I; A0[rz, u] = 1j * Wp * (D2 / a + a * I)   # (P5)
        A0[zz, zz] = -I; A0[zz, u] = -2.0 * Wp * D1       # (P8)
        A0[rt, rt] = -I; A0[rt, v] = Wp * D1              # (P4)
        A0[tt, tt] = -I                                   # (P6)
        A0[tz, tz] = -I; A0[tz, v] = 1j * Wp * a * I      # (P7)
        # ---- premier harmonique (fond : V', T0, T0') ----
        Ac = np.zeros_like(A0); As = np.zeros_like(A0)
        dV_c, dV_s = np.diag(self.V1p), np.diag(self.V2p)
        Sc_, Ss_ = np.diag(self.S1), np.diag(self.S2)
        Scp_, Ssp_ = np.diag(self.S1p), np.diag(self.S2p)
        Ac[rt, rr] = Wp * dV_c;       As[rt, rr] = Wp * dV_s          # (P4)
        Ac[rt, u] = Wp ** 2 * (Sc_ @ D1 - Scp_)
        As[rt, u] = Wp ** 2 * (Ss_ @ D1 - Ssp_)
        Ac[tt, v] = 2.0 * Wp ** 2 * Sc_ @ D1                          # (P6)
        As[tt, v] = 2.0 * Wp ** 2 * Ss_ @ D1
        Ac[tt, rt] = 2.0 * Wp * dV_c; As[tt, rt] = 2.0 * Wp * dV_s
        Ac[tz, u] = 1j * Wp ** 2 / a * Sc_ @ D2                       # (P7)
        As[tz, u] = 1j * Wp ** 2 / a * Ss_ @ D2
        Ac[tz, rz] = Wp * dV_c;       As[tz, rz] = Wp * dV_s
        # ---- terme F' de (P6) : moyenne dans A0, second harmonique ----
        A2c = np.zeros_like(A0); A2s = np.zeros_like(A0)
        A0[tt, u] = A0[tt, u] - Wp ** 3 * np.diag(self.F0p)
        A2c[tt, u] = -Wp ** 3 * np.diag(self.F1p)
        A2s[tt, u] = -Wp ** 3 * np.diag(self.F2p)
        # ---- inertie convective linearisee (campagne G) ----
        if CORIOLIS:
            Ac[u, v] = Ac[u, v] - 2.0 * a ** 2 * C * np.diag(self.V1)   # -2k^2 C V v
            As[u, v] = As[u, v] - 2.0 * a ** 2 * C * np.diag(self.V2)
            g1, g2 = self.V1p, self.V2p                               # -C (V' [+ eps V]) u
            if EPSV:
                g1 = g1 + EPS * self.V1; g2 = g2 + EPS * self.V2
            Ac[v, u] = Ac[v, u] - C * np.diag(g1)
            As[v, u] = As[v, u] - C * np.diag(g2)
        return M, A0, Ac, As, A2c, A2s

    def _reduce(self, M, A):
        Z, keep = self.Z, self.keep
        Mr = (M @ Z)[keep, :]
        Ar = (A @ Z)[keep, :]
        try:
            return solve(Mr, Ar)
        except np.linalg.LinAlgError:
            return lstsq(Mr, Ar, rcond=None)[0]

    def sigma(self, Wp, steady=None):
        """Exposant de Floquet dominant (taux en unites t1) ; steady : valeur propre."""
        st = self.steady if steady is None else steady
        M, A0, Ac, As, A2c, A2s = self._blocks(Wp)
        if st:
            return _sigma_filtered(self._reduce(M, A0 + Ac + A2c), self.Z, self.n)
        if self.De < DE_SW:
            # Symetrie de demi-periode : V, V', T0_rt changent de signe en
            # phi -> phi + pi, F' (second harmonique) non ; l'involution
            # v, T_rt, T_tz -> -v, -T_rt, -T_tz conjugue A(phi+pi) a A(phi)
            # et laisse u (donc le filtre) invariant. Le taux gele est donc
            # pi-periodique : NPH/2 phases de [0, pi) donnent exactement la
            # moyenne de NPH phases de [0, 2pi) (verifie a 1e-15, 16/09).
            acc = 0.0
            phs = np.linspace(0.0, np.pi, NPH // 2, endpoint=False)
            for ph in phs:
                B = self._reduce(M, A0 + np.cos(ph) * Ac + np.sin(ph) * As
                                 + np.cos(2 * ph) * A2c + np.sin(2 * ph) * A2s)
                acc += _sigma_filtered(B, self.Z, self.n)
            return acc / len(phs)
        w1 = self.De / Wp
        T1 = 2.0 * np.pi / w1
        dt = T1 / self.m
        Phi = np.eye(8 * self.n - 6, dtype=complex)
        logs = 0.0
        for j in range(self.m):
            ph = w1 * (j + 0.5) * dt
            B = self._reduce(M, A0 + np.cos(ph) * Ac + np.sin(ph) * As
                             + np.cos(2 * ph) * A2c + np.sin(2 * ph) * A2s)
            Phi = expm(B * dt) @ Phi
            nrm = np.max(np.abs(Phi))
            if np.isfinite(nrm) and nrm > 0:
                Phi /= nrm
                logs += np.log(nrm)
        if not np.all(np.isfinite(Phi)):
            return 1e3
        try:
            mu, V = eig(Phi)
        except np.linalg.LinAlgError:
            return 1e3
        for i in np.argsort(np.abs(mu))[::-1]:
            if not np.isfinite(mu[i]) or abs(mu[i]) == 0.0:
                continue
            if _tail_fraction(self.Z @ V[:, i], self.n) < TAIL:
                return float((np.log(abs(mu[i])) + logs) / T1)
        return -1e3


# ====================================================================== balayage
DG = 0.01
G0, G1 = 0.10, 7.75
DA = 0.5
A_LO, A_HI = 2.0, 40.0
HALF0 = 1.0
HALF_MAX = 6.0
FULL_EVERY = 50
RTOL = 3e-3


def sb_of(S):
    """S = 1 est regulier pour l'operateur complet : Sb = 0 exactement."""
    return 0.0 if S >= 1.0 else (1.0 - S) / S


def wp_hi_of(E):
    """Plafond du balayage log en Wp. Releve (16/09) : a gamma = 20 la campagne
    F' atteint Wp = 1.4e3 (E=1), 1.4e4 (E=10), 6.9e4 (E=50) pour S = 0.3."""
    if E >= 50.0:
        return 5.0e5
    if E >= 10.0:
        return 1.0e5
    if E >= 1.0:
        return 1.0e4
    return WP_HI


def thresh(fb, hint=None, hi=None):
    """Seuil de fb.sigma ; bracket serre autour de hint, repli sur scan log."""
    g = fb.sigma
    hi = WP_HI if hi is None else hi
    if hint is not None and np.isfinite(hint):
        for fac in (1.06, 1.3, 2.2):
            a, b = max(WP_LO, hint / fac), min(hi, hint * fac)
            if a >= b:
                continue
            sa = g(a)
            if sa >= 0:
                continue
            sb_ = g(b)
            if sb_ >= 0:
                try:
                    return brentq(g, a, b, xtol=1e-4, rtol=RTOL)
                except Exception:
                    return 0.5 * (a + b)
    grid = np.geomspace(WP_LO, hi, NSCAN)
    pw, ps = grid[0], g(grid[0])
    if ps >= 0:
        return np.nan
    for W in grid[1:]:
        s = g(W)
        if ps < 0 <= s:
            try:
                return brentq(g, pw, W, xtol=1e-4, rtol=RTOL)
            except Exception:
                return 0.5 * (pw + W)
        pw, ps = W, s
    return np.nan


def coarse_seed(S, E, gam):
    """(k_c, Wp_c) d'amorcage : SEED injecte, sinon campagne F' locale."""
    if SEED:
        k = "%.3f" % gam
        if k in SEED:
            return float(SEED[k][0]), float(SEED[k][1])
        ks = sorted(SEED, key=lambda s: abs(float(s) - gam))
        return float(SEED[ks[0]][0]), float(SEED[ks[0]][1])
    for f in (os.path.join("resultats_F", "fine_S%.2f_E%g.csv" % (S, E)),
              os.path.join("resultats_F", "fine_S%.2f_E%g.csv" % (0.98 if S >= 1 else S, E))):
        if os.path.exists(f):
            d = np.atleast_2d(np.genfromtxt(f, delimiter=",", skip_header=1))
            d = d[np.isfinite(d[:, 1])]
            if len(d):
                i = int(np.argmin(np.abs(d[:, 0] - gam)))
                return float(d[i, 2]), float(d[i, 1])
    return None, None


def scan(Sb, E, gam, alist, prev, hi=None):
    out = {}
    hint = None
    for a in alist:
        h = prev.get(a, hint)
        w = thresh(ShaqfehFull(a, Sb, E, gam), hint=h, hi=hi)
        if np.isfinite(w):
            out[a] = w
            hint = w
    return out


def best_of(d):
    if not d:
        return np.nan, np.nan
    a = min(d, key=d.get)
    return d[a], a


NAN_RETRY = 10          # balayage complet sans seuil : ne le relancer qu'apres NAN_RETRY points


def march(S, E, g0, g1, out, gams=None):
    """Marche en gamma, reprenable (saute les gamma deja ecrits dans out).
    gams : liste explicite (basse frequence, queues) ; sinon arange(g0, g1, DG).
    Chaque point est aussi imprime (ligne ROW) : un kernel coupe a 12 h ne
    committe que son journal, d'ou l'on peut alors recuperer les donnees."""
    Sb = sb_of(S)
    HI = wp_hi_of(E)
    if gams is None:
        gams = np.round(np.arange(g0, g1 + 1e-9, DG), 3)
    else:
        gams = np.round(np.asarray(gams, dtype=float), 5)
    done = set()
    for src in [out] + glob.glob("fine_S%.2f_E%g_pfix.csv" % (S, E)):
        if not os.path.exists(src):
            continue
        for line in open(src).read().splitlines()[1:]:
            if line.strip():
                done.add(round(float(line.split(",")[0]), 5))
    fh = open(out, "a")
    if not done:
        fh.write("gamma, Wp_c, alpha_c, Wi_c_eq\n"); fh.flush()

    def emit(gam, w, k):
        if np.isfinite(w):
            line = "%.5f, %.6g, %g, %.6g" % (gam, w, k, w / np.sqrt(EPS))
        else:
            line = "%.5f, nan, nan, nan" % gam
        fh.write(line + "\n"); fh.flush()
        print("ROW S=%g E=%g %s" % (S, E, line), flush=True)

    kc, prev = None, {}
    kc, w0 = coarse_seed(S, E, float(gams[0]))
    if kc is not None and np.isfinite(w0):
        prev = {kc: w0}
    t0, nrun, last_nan = time.time(), 0, None
    for n, gam in enumerate(gams):
        if round(float(gam), 5) in done:
            continue
        if kc is None or (n % FULL_EVERY == 0 and nrun > 0):
            if kc is None and last_nan is not None and n - last_nan < NAN_RETRY:
                emit(gam, np.nan, np.nan)      # aucun mode resolu tout pres : on n'insiste pas
                nrun += 1
                continue
            coarse = np.round(np.arange(A_LO, A_HI + 1e-9, 2.0), 2)
            dc = scan(Sb, E, gam, coarse, prev, hi=HI)
            _, kk = best_of(dc)
            if not np.isfinite(kk):
                emit(gam, np.nan, np.nan)
                kc = None; last_nan = n; nrun += 1
                continue
            kc = kk
            prev = dc
            last_nan = None
        half = HALF0
        d = {}
        wmin, kmin = np.nan, np.nan
        while True:
            alist = np.round(np.arange(max(A_LO, kc - half),
                                       min(A_HI, kc + half) + 1e-9, DA), 2)
            d = scan(Sb, E, gam, alist, prev, hi=HI)
            wmin, kmin = best_of(d)
            if not np.isfinite(wmin):
                break
            at_edge = (kmin <= alist[0] + 1e-9 and alist[0] > A_LO) or \
                      (kmin >= alist[-1] - 1e-9 and alist[-1] < A_HI)
            if kmin >= A_HI - 1e-9:
                sys.stderr.write("ATTENTION k_c au plafond A_HI=%g a gamma=%.3f "
                                 "(S=%g E=%g) : seuil = MAJORANT\n" % (A_HI, gam, S, E))
            if not at_edge or half >= HALF_MAX:
                break
            half += 2.0
            kc = kmin
        prev = d if d else prev
        if np.isfinite(wmin):
            kc = kmin
            emit(gam, wmin, kmin)
        else:
            kc = None
            emit(gam, np.nan, np.nan)
        nrun += 1
    fh.close()
    return time.time() - t0, nrun


def task(t):
    S, E, g0, g1, k = t
    out = "fine_S%.2f_E%g_p%d.csv" % (S, E, k)
    dt, n = march(S, E, g0, g1, out)
    return ("S=%-4g E=%-5g bloc %d : %d gammas en %.0f s (%.1f s/gamma) -> %s"
            % (S, E, k, n, dt, dt / max(n, 1), out))


def _selftest():
    """(1) egalite avec le prototype du debat ; (2) Taylor stationnaire ;
    (3) S = 1 regulier ; (4) chronometrage."""
    t0 = time.time()
    print("== self-test run_full_sweep (operateur complet) ==")
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from _coriolis_test import ShaqfehWFC
        for (S, E, g, k, wp) in ((0.98, 1.0, 1.256, 6.0, 1.7866), (0.3, 10.0, 0.45, 6.0, 394.21)):
            Sb = (1 - S) / S
            s1 = ShaqfehFull(k, Sb, E, g).sigma(wp)
            s2 = ShaqfehWFC(k, Sb, E, g, coriolis=True).sigma(wp)
            print("  (1) S=%.2f E=%g g=%.3f k=%g Wp=%.4g : sigma=%+.6e  prototype=%+.6e  diff=%.1e"
                  % (S, E, g, k, wp, s1, s2, abs(s1 - s2)), flush=True)
    except Exception as e:
        print("  (1) prototype _coriolis_test indisponible :", e)
    fb = ShaqfehFull(3.12, 100.0, 0.02435, 1.0, steady=True)
    print("  (2) Taylor stationnaire (Sb=100, k=3.12, E=0.02435 -> Ta=3441) : sigma(Wp=1)=%+.4f (attendu ~0)"
          % fb.sigma(1.0), flush=True)
    g, E, k = 1.256, 1.0, 6.0
    t1 = time.time()
    s = ShaqfehFull(k, 0.0, E, g).sigma(1.8)
    print("  (3) S=1 (Sb=0) : sigma(Wp=1.8)=%+.4f  [%.1f s / evaluation, N=%d]" % (s, time.time() - t1, N), flush=True)
    w = thresh(ShaqfehFull(k, 0.0, E, g), hint=1.8)
    print("      seuil S=1 E=1 gamma=1.256 k=6 : Wp_c=%.4g  Wi_c=%.4g" % (w, w / np.sqrt(EPS)), flush=True)
    t1 = time.time()
    ShaqfehFull(7.0, 0.0, 0.01, 10.0, n=48).sigma(2.0)
    print("  (4) N=48 (queue) : %.1f s / evaluation" % (time.time() - t1))
    print("self-test %.0f s" % (time.time() - t0))


def main():
    nums = []
    for a in sys.argv[1:]:
        try:
            nums.append(float(a))
        except ValueError:
            pass
    if "--test" in sys.argv:
        _selftest()
        return
    if len(nums) >= 5:
        print(task((nums[0], nums[1], nums[2], nums[3], int(nums[4]))))
        return
    print(__doc__)


if __name__ == "__main__":
    main()

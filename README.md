# Oscillatory Taylor–Couette flow of an Oldroyd-B fluid — stability codes

Linear stability of the **modulated (co-oscillating) Taylor–Couette flow** of an
Oldroyd-B fluid in the narrow-gap limit. Both cylinders oscillate azimuthally
about a zero mean, Ω_in = Ω₁cos ωt and Ω_out = Ω₂cos ωt. The codes give the
critical Weissenberg number Wi_c and the critical axial wavenumber k_c as
functions of the dimensionless frequency γ, for viscosity ratios S = η_p/η and
elasticity numbers E = We/Re.

Only the codes needed to reproduce the calculation are kept here.

## Python — current formulation

The narrow-gap ordering is taken at fixed elasticity number. The inertia of the
base flow and of the disturbance is kept at leading order in the gap-to-radius
ratio ε, and the problem closes on **eight fields**: the radial and azimuthal
velocities and the six polymer stresses. The azimuthal momentum balance carries
the transport (∂ₓV⁰ + 2εV⁰)u, that is the full Coriolis pair of the
co-oscillating base. It matters below De = 0.3 and is included there.

### Solver

| File | Role |
|---|---|
| `python/run_full_sweep.py` | the eight-field operator, the Floquet monodromy, the low-frequency cycle-average branch, and the marching sweep in γ with a moving window in k |

```bash
python run_full_sweep.py 0.9 1 0.10 2.00 1   # S, E, gamma range, block number
```

### Assembling the curves

| File | Role |
|---|---|
| `python/_G_quicklook.py` | merges the blocks of a sweep into one CSV per case, `G_S<S>_E<E>.csv` |
| `python/_H_merge.py` | same for the low-frequency recomputation with the full Coriolis pair, `H_S<S>_E<E>.csv` |
| `python/_GH_merge.py` | assembles the curve used in the article: H below De = 0.3, G above, `GH_S<S>_E<E>.csv` |
| `python/_kscan_merge.py` | where two branches compete, replaces the windowed minimum by the true minimum of a full scan in k (k = 1 … 80). A minimum on the edge of the grid is refused |

```bash
python _G_quicklook.py
python _H_merge.py
python _GH_merge.py
python _kscan_merge.py
```

### Tables and figures

| File | Role |
|---|---|
| `python/_GH_table.py` | critical parameters of each case: γ*, De*, Wi_c*, k_c*, Wi_c(0.2), Wi_c(20) |
| `python/_gen_wkb_table.py` | constants A(S) of the low-frequency law Wi_c → ε^(−1/2) A(S) γ^(−2) of the inertialess system, read from the quasi-static campaign |
| `python/_G_figs_article.py` | one figure per case: Wi_c(γ) in black on the left axis, k_c(γ) in blue on the right axis, the inertialess law dashed in red, LaTeX rendering |
| `python/_G_inserer_figures.py` | puts the figures into the LaTeX source |

```bash
python _GH_table.py
python _gen_wkb_table.py
python _G_figs_article.py
python _G_inserer_figures.py
```

What the figure script does, and says in the caption:

* the γ axis is graduated 0.1, 0.5, 1, 5, 10, 20;
* at E = 0.01 the computation stops at γ = 0.18, and both curves are continued
  down to γ = 0.1 by their local power law, anchored on the first computed
  point (value and slope continuous);
* over 10 ≤ γ ≤ 20 the plotted curves are smoothed for display only, by a
  running median over 11 points followed by a 5-point mean. The data files and
  the tables are the raw values;
* k_c is dashed where it sits on the edge of the wavenumber grid: the minimum
  over k is not attained there and Wi_c is only an upper bound.

## MATLAB — Floquet codes

| File | Role |
|---|---|
| `matlab/sigma_B_floquet.m` | largest Floquet exponent of the modulated system at given Wi, k, frequency, E and S |
| `matlab/B_freq_for_one_S.m` | sweep in frequency at one viscosity ratio |
| `matlab/B_for_one_S.m` | critical values Wi_c(γ), k_c(γ) at one viscosity ratio |
| `matlab/B_loop_over_S.m` | loop over the viscosity ratios |
| `matlab/chebdif.m` | Chebyshev differentiation matrices (Weideman & Reddy) |

```matlab
B_for_one_S(0.9)
```

## Numerics

Chebyshev collocation on the Gauss–Lobatto points, N = 28 over most of the
frequency range and N = 48 in the high-frequency tail. The forcing period is
cut into m = 40 equal steps; on each step the operator is frozen at the
midpoint and advanced by its matrix exponential; the product of the steps is
the monodromy matrix, renormalised at every step. A mode is discarded when more
than 5 % of the energy of its radial velocity lies in the upper third of its
Chebyshev spectrum. The threshold at fixed k is a root of the Floquet exponent,
found by Brent's method. The critical values are the minimum over a window in k
that follows k_c from one frequency to the next, and the whole range of k is
scanned again every 50 frequencies. Below De = 0.05 the period is too long for
the matrix exponentials, and the threshold comes from the cycle average of the
largest frozen growth rate, taken over 24 phases of [0, π) by the half-period
symmetry of the base flow.

Resolution check at low frequency: at E = 0.01 and γ ≤ 0.14, N = 28 is not
converged (at S = 0.7, k = 57, γ = 0.14 the threshold is 3.7·10³ at N = 28 and
1.35·10² at N = 56; N = 56 and N = 68 agree to all printed digits). That is why
the E = 0.01 curves are computed from γ = 0.18.

Data files are not included. The plotting script reads one CSV per case,
`resultats_G/GH_S<S>_E<E>.csv`, with the columns `gamma, Wp_c, k_c, Wi_c`,
where Wp = ε^(1/2) Wi and ε = 0.14.

## Known limitations

* Where k_c reaches the edge of the wavenumber grid (k = 40 over most of the
  range, k = 80 in the tail and in the full scans), the minimum over k is not
  attained and the threshold is only an upper bound. This happens on the
  high-frequency branch of the near-UCM case S = 0.98 at E ≥ 1.
* Below γ ≈ 0.15 the narrow-gap model itself is at its limit (γ² comparable
  to ε²), and N = 28 is not converged there.
* The cycle-average branch used below De = 0.05 is the exponent of a
  continuously followed branch only while the same branch stays the most
  unstable over the whole cycle.

## References

M. Hayani Choujaa, M. Riahi, S. Aniss, *Floquet-based analysis on
three-dimensional non-axisymmetric instabilities in oscillatory-driven
Taylor–Couette flows and their low-frequency asymptotic behavior using
Wentzel–Kramers–Brillouin method*, Phys. Fluids **36**, 014101 (2024).

M. Hayani Choujaa, S. Aniss, M. Ouazzani Touhami, *Stability of an oscillatory
Taylor–Couette flow in an upper-convected Maxwell fluid*, Phys. Fluids **33**,
074105 (2021).

The Oldroyd-B study that these codes serve is in preparation.

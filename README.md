# Oscillatory Taylor–Couette flow of an Oldroyd-B fluid — stability codes

Linear stability of the **modulated (co-oscillating) Taylor–Couette flow** of an
Oldroyd-B fluid in the narrow-gap limit. Both cylinders oscillate azimuthally
about a zero mean, Ω_in = Ω₁cos ωt and Ω_out = Ω₂cos ωt, and the codes compute
the critical Weissenberg number Wi_c and the critical axial wavenumber k_c as
functions of the dimensionless frequency γ, for several viscosity ratios
S = η_p/η and elasticity numbers E = We/Re.

Two generations of code live here.

* **MATLAB** — the earlier Floquet codes (Chebyshev collocation, monodromy over
  one forcing period), used for the UCM and Oldroyd-B systems.
* **Python** — the current formulation: the narrow-gap ordering is taken at
  fixed elasticity number, so the inertia of the base flow and of the
  disturbance is kept at leading order in the gap-to-radius ratio, and the
  problem closes on **eight fields** (two velocity components and the six
  polymer stresses). Floquet monodromy at moderate and high frequency, WKB
  cycle average at low frequency.

## Layout

| Path | Contents |
|---|---|
| `matlab/floquet/` | `sigma_B_floquet.m` (largest Floquet exponent), `B_for_one_S.m`, `B_freq_for_one_S.m`, `B_loop_over_S.m` (sweeps in γ, k, S), `characterize_modes.m`, `replot_modes.m`, `plot_S_09_fine.m`, `chebdif.m` (Weideman–Reddy differentiation matrices) |
| `matlab/purely_elastic/` | the purely elastic runs at fixed E (E = 0.006, 0.01, 1, 10) and `run_S07_matlab.m` |
| `python/run_full_sweep.py` | the eight-field operator, the monodromy, the WKB branch, the marching sweep in γ with the window in k |
| `python/_gen_G_kaggle.py` | splits a full sweep (γ = 0.01 → 20, 16 cases) into independent blocks for a batch service |
| `python/_gen_G_complement.py` | regenerates only the frequencies missing from a sweep, in short blocks |
| `python/_G_recup_logs.py` | recovers computed points from the logs of blocks that were cut before they finished |
| `python/_G_quicklook.py` | merges the blocks into one CSV per case and draws a quick comparison figure |
| `python/_G_figs_article.py` | the article figures: one figure per case, Wi_c(γ) in black on the left axis, k_c(γ) in blue on the right axis, LaTeX rendering |
| `python/_G_apercu.py` | one overview figure, four elasticity numbers × four viscosity ratios |
| `python/_G_inserer_figures.py` | puts a figure in the LaTeX manuscript only when its case is computed over the whole frequency range, and an empty frame otherwise |

## Running

```bash
# one block of a sweep: S, E, gamma range, block number
python run_full_sweep.py 0.9 1 0.10 2.00 1

# merge the blocks, then draw
python _G_quicklook.py
python _G_figs_article.py
```

```matlab
% one viscosity ratio, sweep in frequency
B_for_one_S(0.9)
```

The Python scripts read one CSV per case, `resultats_G/G_S<S>_E<E>.csv`, with
the columns `gamma, Wp_c, k_c, Wi_c`. The data files themselves are not in this
repository.

## Numerics

Chebyshev collocation (N = 28 in the body of the sweep, N = 48 in the
high-frequency tail). The period is cut into m = 40 equal steps; on each step
the operator is frozen at the midpoint and advanced by its matrix exponential,
and the product of the steps is the monodromy matrix, renormalised at every
step. A mode is discarded when more than 5 % of the energy of its radial
velocity lies in the upper third of its Chebyshev spectrum. The threshold at
fixed k is a root of the Floquet exponent, found by Brent's method, and the
critical values are the minimum over a window in k that follows k_c from one
frequency to the next. Below De = 0.05 the period is too long for the matrix
exponentials, and the threshold comes instead from the cycle average of the
largest frozen growth rate, taken over 24 phases of [0, π) by the half-period
symmetry of the base flow.

## Known limitations

* The narrow-gap truncation keeps the radial half of the Coriolis coupling of
  the co-oscillating base and drops the azimuthal half. Measured effect on the
  threshold: at most 2 % for De ≥ 0.3, but 20 % or more for De ≤ 0.09. The
  low-frequency branch is therefore under revision.
* Where k_c reaches the edge of the wavenumber grid (k = 40 in the body of the
  sweep, k = 80 in the tail), the minimum over k is not attained and the
  threshold is only an upper bound.
* The WKB branch uses the cycle average of the largest frozen growth rate. That
  average equals the exponent of a continuously followed branch only while the
  same branch stays the most unstable over the whole cycle.

## References

M. Hayani Choujaa, M. Riahi, S. Aniss, *Floquet-based analysis on
three-dimensional non-axisymmetric instabilities in oscillatory-driven
Taylor–Couette flows and their low-frequency asymptotic behavior using
Wentzel–Kramers–Brillouin method*, Phys. Fluids **36**, 014101 (2024).

M. Hayani Choujaa, S. Aniss, M. Ouazzani Touhami, *Stability of an oscillatory
Taylor–Couette flow in an upper-convected Maxwell fluid*, Phys. Fluids **33**,
074105 (2021).

The Oldroyd-B study that these codes serve is in preparation.

# Oscillatory Taylor–Couette flow of an Oldroyd-B fluid — stability codes

Linear stability of the **modulated (co-oscillating) Taylor–Couette flow** of an
Oldroyd-B fluid in the narrow-gap limit. Both cylinders oscillate azimuthally
about a zero mean, Ω_in = Ω₁cos ωt and Ω_out = Ω₂cos ωt. The codes give the
critical Weissenberg number Wi_c and the critical axial wavenumber k_c as
functions of the dimensionless frequency γ, for viscosity ratios S = η_p/η and
elasticity numbers E = We/Re.

Only the codes needed to reproduce that calculation are kept here.

## Python — current formulation

The narrow-gap ordering is taken at fixed elasticity number, so the inertia of
the base flow and of the disturbance is kept at leading order in the
gap-to-radius ratio, and the problem closes on **eight fields**: the radial and
azimuthal velocities and the six polymer stresses.

| File | Role |
|---|---|
| `python/run_full_sweep.py` | the eight-field operator, the Floquet monodromy, the low-frequency WKB branch, and the marching sweep in γ with a moving window in k |
| `python/_G_quicklook.py` | merges the blocks of a sweep into one CSV per case, `G_S<S>_E<E>.csv` |
| `python/_G_figs_article.py` | one figure per case: Wi_c(γ) in black on the left axis, k_c(γ) in blue on the right axis, LaTeX rendering |

```bash
python run_full_sweep.py 0.9 1 0.10 2.00 1   # S, E, gamma range, block number
python _G_quicklook.py                       # merge the blocks
python _G_figs_article.py                    # draw
```

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

Chebyshev collocation, N = 28 over most of the frequency range and N = 48 in the
high-frequency tail. The forcing period is cut into m = 40 equal steps; on each
step the operator is frozen at the midpoint and advanced by its matrix
exponential, and the product of the steps is the monodromy matrix, renormalised
at every step. A mode is discarded when more than 5 % of the energy of its
radial velocity lies in the upper third of its Chebyshev spectrum. The threshold
at fixed k is a root of the Floquet exponent, found by Brent's method, and the
critical values are the minimum over a window in k that follows k_c from one
frequency to the next. Below De = 0.05 the period is too long for the matrix
exponentials, and the threshold comes from the cycle average of the largest
frozen growth rate, taken over 24 phases of [0, π) by the half-period symmetry
of the base flow.

Data files are not included. The plotting script reads one CSV per case,
`resultats_G/G_S<S>_E<E>.csv`, with the columns `gamma, Wp_c, k_c, Wi_c`.

## Known limitations

* The narrow-gap truncation keeps the radial half of the Coriolis coupling of
  the co-oscillating base and drops the azimuthal half. Measured effect on the
  threshold: at most 2 % for De ≥ 0.3, but 20 % or more for De ≤ 0.09, so the
  low-frequency branch is under revision.
* Where k_c reaches the edge of the wavenumber grid (k = 40 over most of the
  range, k = 80 in the tail), the minimum over k is not attained and the
  threshold is only an upper bound.
* The WKB branch uses the cycle average of the largest frozen growth rate. That
  average is the exponent of a continuously followed branch only while the same
  branch stays the most unstable over the whole cycle.

## References

M. Hayani Choujaa, M. Riahi, S. Aniss, *Floquet-based analysis on
three-dimensional non-axisymmetric instabilities in oscillatory-driven
Taylor–Couette flows and their low-frequency asymptotic behavior using
Wentzel–Kramers–Brillouin method*, Phys. Fluids **36**, 014101 (2024).

M. Hayani Choujaa, S. Aniss, M. Ouazzani Touhami, *Stability of an oscillatory
Taylor–Couette flow in an upper-convected Maxwell fluid*, Phys. Fluids **33**,
074105 (2021).

The Oldroyd-B study that these codes serve is in preparation.

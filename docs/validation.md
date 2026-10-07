# Verification protocol

All experiments are dimensionless. The two materials are invented; they are
not measured material properties. An initial temperature excess of zero is
the baseline, so negative values mean an unphysical undershoot below that
baseline, not a claim about negative absolute temperature.

1. **Independent algebra:** compare theta steps to a dense solve, spectral
   propagation to SciPy's dense matrix exponential, and all three methods to
   the closed-form two-cell exchange mode. Check exact paired face resistance,
   symmetry, heat preservation, constant preservation, and the CN energy identity.
2. **Pulse stress test:** 128 cells on [0,1], conductivity 1/4 and volumetric
   capacity 1/2 on the left/right halves; a unit pulse in zero-based cell 63.
   Use dt = 20 times min(M_i/K_ii), 12 macro steps. Rannacher uses four BE
   half-steps over the first two macro intervals. Record all 13 states for
   each method. Tests assert decreasing quadratic energy, small heat drift,
   the CN undershoot and the observed Rannacher error reduction.
3. **Time error:** same fixed mesh and initial pulse, final time 0.02, with
   20, 40, 80, 160 and 320 macro intervals. Compare with exp(-M^-1 K t)T0 in
   the mass-weighted RMS norm. This is a discrete-ODE reference; it says
   nothing about continuum accuracy at the material interface.
4. **Space error:** homogeneous k=c=1, initial field 2+cos(pi*x), final time
   0.1. Initialize exact cell averages, propagate the discrete equations
   spectrally, and compare to cell averages of 2+exp(-pi^2*t)cos(pi*x).
   Use 16,32,64,128,256 uniform cells. This measures continuum cell-average
   error, not reconstructed-field error. It does not establish interface
   convergence for heterogeneous or arbitrary nonuniform meshes.

The final pair of refinement levels defines each reported order. Full error
sequences are retained because coarse-step CN errors have not reached their
second-order regime. Rannacher is not more accurate than CN at every step size.

`python scripts/check_samples.py` reproduces every JSON number and every CSV
value in a temporary directory. Metric tolerances are rtol=5e-5, atol=5e-12
(including near-roundoff drift); history tolerances are rtol=2e-10,
atol=5e-12. Figures are regenerated but not tested by byte equality across
Matplotlib versions. Manuscript inputs and its delivered PDF have SHA-256
hashes in `paper/validation.json`; `scripts/check_manifest.py` verifies them.

CI runs tests, Ruff, full numerical reproduction and wheel/source builds on
Ubuntu/Python 3.11 and Windows/Python 3.14, then tests the installed wheel.
A third job compiles the manuscript with pdfLaTeX, rejects undefined references
and overfull boxes, checks the five-page limit and renders every page for review.
No timing benchmark or performance advantage is claimed.

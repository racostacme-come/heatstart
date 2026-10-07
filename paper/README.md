# Building and reviewing the manuscript

From the repository root, install `.[dev]` and either Tectonic or TeX Live.

```sh
python scripts/check_samples.py
python paper/build.py --engine tectonic
# Alternative: python paper/build.py --engine pdflatex
python scripts/check_manifest.py
pdftoppm -r 110 -png paper/paper.pdf paper/build/page
```

The build script derives `numbers.tex` from recorded metrics, compiles the
LaTeX (two passes for pdfLaTeX), rejects undefined references and overfull
boxes, checks all pages are nonempty, and enforces a maximum of five pages
including references. The expected manuscript has four pages. It writes the
deliverable `paper.pdf` and a hash manifest `validation.json`. Figures come
only from this repository's synthetic computations.

Numerical reproduction runs independently before compilation; generating
LaTeX macros alone does not validate the numerical results. After any input
change, regenerate results if needed, rebuild, run the manifest check and
visually inspect every rendered page. The manuscript CI artifact includes
the independently compiled PDF, page PNGs, log and manifest for review.

The four-page length stated in the README describes this version, whereas
the automated acceptance gate enforces the series requirement of at most
five pages. Private course files are not required to reproduce the paper.

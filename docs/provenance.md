# Concept provenance and originality

On 7 October 2026 the three DTU source folders (ATSA, MLSP, MMME) were inspected
read-only. The selected source was MMME course 41747, **Modelling in Materials
and Manufacturing Engineering**. The actual Lecture 3 notes discuss heat
diffusion, material interfaces and analytical thermal contact; Lecture 4 covers
explicit Fourier-number stability, positive update weights and series half-cell
thermal resistances. Those concepts motivate this project's heterogeneous
finite-volume model and its separate positivity and energy diagnostics.

This project extends those concepts to implicit time stepping and startup
damping. It is an original implementation and synthetic experiment, not a
transcription or solution of a course exercise. No lecture text, figures,
exercise implementations, exams or private data are distributed. No DTU or
ATSA files were modified. The course is not an affiliation or endorsement.

The current inspection also identified new ATSA Lecture 6 material on
hierarchical forecast reconciliation (including heat-load forecasting) and
probabilistic forecast evaluation (PIT, scoring rules and multivariate
evaluation), added on 7 October. MLSP chapter 5 discusses sampling and
quantization; this was the primary subject of the previous AliasGuard project.
These inspections support rotation without guessing topics from acronyms.

The series already uses C++17 and Python in several projects. HeatStart uses
Python with SciPy's compiled banded linear algebra: the original scientific
work is the discretization, error separation and verification, and a second
custom C++ solver would duplicate the existing linear algebra library.

## Public references

- R. Rannacher, *Finite element solution of diffusion problems with irregular
  data*, Numerische Mathematik **43**, 309–327 (1984).
  [Bibliographic record](https://eudml.org/doc/132904).
- H. P. Langtangen, *Finite difference methods for diffusion processes*,
  INF5620 teaching notes (2014), especially amplification-factor analysis and
  startup damping.
  [Author's public notes](https://hplgit.github.io/INF5620/doc/pub/H14/diffu/html/._main_diffu001.html).

References describe established methods. No algorithmic novelty, external
review, experimental validation, or institutional authorship is claimed.

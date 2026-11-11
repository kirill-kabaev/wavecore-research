# Publications

This page lists WaveCore Research manuscripts and technical research documents.
The entries below are repository documents; they are not presented as
peer-reviewed journal publications.

## Physical Foundations of Numerical Modeling of Electromagnetic Waves

**Author:** Kirill Kabaev<br>
**Year:** 2026 (the manuscript uses `\\today`)<br>
**Status:** Foundational research manuscript / not claimed as peer reviewed

This manuscript introduces the physical and mathematical background for
computational electromagnetics. It reviews Maxwell's equations, wave and
Helmholtz formulations, time-harmonic conventions, and boundary conditions,
then connects those models with numerical discretization and sparse linear
algebra. The later sections discuss numerical methods, sparse solvers, and
telecommunications-oriented applications as research context. It is a
foundational draft for the WaveCore project rather than a report of a completed
computational validation campaign. The document is written in Russian and is
maintained as a modular LuaLaTeX source with BibLaTeX references.

**Research themes:** Maxwell equations; Helmholtz equations; computational
electromagnetics; boundary conditions; numerical methods; sparse solvers.

- LaTeX source: <a href="https://github.com/kirill-kabaev/wavecore-research/blob/58fcf252f97e69b60b6bf1efe2fd0ece9286c8a6/papers/01_maxwell_helmholtz/main.tex">`papers/01_maxwell_helmholtz/main.tex`</a>
- [Locally built PDF](papers/maxwell-helmholtz.pdf)

## Numerical Robustness of Restarted GMRES for Complex-Valued Linear Systems: Preliminary Computational Experiments

**Author:** Kirill Kabaev, WaveCore Research<br>
**Year:** 2026<br>
**Status:** Preliminary technical research note / Not peer reviewed

This preliminary technical note reports reproducible numerical checks of a
sequential, unpreconditioned restarted GMRES implementation for complex sparse
linear systems. It describes Arnoldi iteration with the Hermitian inner
product, two-pass modified Gram--Schmidt orthogonalization, complex Givens
rotations, and an explicit true-residual stopping check. Deterministic tests
probe nonnormal and ill-conditioned systems, restart behavior, stagnation,
breakdown, finite scaling, Arnoldi orthogonality, and disagreement between
least-squares residual estimates and recomputed residuals. The discussion keeps
forward solution error separate from residual accuracy and records numerical
limitations. It is verification evidence for one implementation, not a claim
of universal stability or scientific novelty.

**Research themes:** restarted GMRES; complex Krylov methods; Arnoldi
orthogonalization; numerical robustness; residual diagnostics; sparse linear
systems.

- LaTeX source: <a href="https://github.com/kirill-kabaev/wavecore-research/blob/58fcf252f97e69b60b6bf1efe2fd0ece9286c8a6/papers/02_gmres_numerical_observations/main.tex">`papers/02_gmres_numerical_observations/main.tex`</a>
- PDF: not currently published on the website; build it locally from the
  LaTeX source with `latexmk -lualatex`.

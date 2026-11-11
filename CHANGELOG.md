# Changelog

This file records verified project checkpoints. Research milestones are not
automatically published releases or peer-reviewed results.

## Unreleased — Stage 5C: verification hardening

- Stage 5C.1, committed as b99f45ae8e056efe2582ddf6186db48b3bf4aafe:
  shared test-only convergence acceptance checks and persistent regression
  fixtures, including typed failure reasons, exact diagnostic conditions,
  floating-point boundaries, and a real iteration-limited GMRES solve.
  This tests acceptance semantics without changing Stage 5B numerical criteria.
- Stage 5C.2, committed as 58fcf252f97e69b60b6bf1efe2fd0ece9286c8a6:
  a separate rectangular convergence CTest and CSV, with independent references,
  weighted norms, spectral-gap checks, true residuals and named finite metrics.
  The suite at this checkpoint contains thirteen CTest registrations.
- Stage 5C.3 documentation: the
  [Numerical Verification Protocol](docs/numerical-verification.md) consolidates
  the Stage 5A/5B/5C.1/5C.2 model, thresholds, references, commands and limitations.
  Historical CSVs, numerical algorithms, and acceptance criteria are preserved.

These are research checkpoints, not published releases. The rectangular study
covers one mode and fixed spacing ratio on one platform; it adds no general
stability, near-resonance, performance or cross-platform guarantee.

## Unreleased — Stage 5B: spatial convergence

### Added

- Separate CTest experiments `wavecore.convergence` and
  `wavecore.convergence.multimode`, using continuous manufactured forcing on
  7x7, 15x15, 31x31, and 63x63 interior-node grids on the unit square, k=1.
- A single sine mode (1,1) and the independent mixture
  phi_11 - 0.3 phi_23, both with zero Dirichlet boundaries.
- Preserved full-precision tables:
  [single-mode CSV](experiments/helmholtz_2d/results/stage5b_grid_convergence.csv)
  and [multimode CSV](experiments/helmholtz_2d/results/stage5b_multimode_convergence.csv).
- Exact discrete modal references, independent true-residual checks, weighted
  L2 and Linf errors, observed orders, and separate numerical/algebraic errors.
- Failure diagnostics containing grid size, spacing, condition name, measured
  quantities, and allowed limits; numerical criteria and CSV results are preserved.
- [Technical decisions](docs/technical-decisions.md), expanded test coverage,
  research-log records, and roadmap status.

### Verified scope

- Both norms exhibit second-order spatial convergence for the recorded cases.
  On the final refinement the orders are approximately 2.000568 / 2.000568
  (single mode) and 2.003139 / 2.002857 (mixture), in L2 / Linf.
- Success requires CONVERGED, a separately recomputed true relative residual
  <=1e-11, and numerical/algebraic errors <=0.1% of spatial errors.
- The mixture takes 2,2,2,3 GMRES iterations in the recorded environment,
  avoiding an exclusively one-step eigenvector verification.
- At the Stage 5B checkpoint the suite contained eleven CTest registrations. Numerical-core and manuscript
  sources are unchanged by Stage 5B.

### Limitations

Four square isotropic grids, two fixed manufactured fields, subresonant k=1,
and the recorded GCC/WSL2 environment. No general stability, near-resonance,
anisotropic-grid, scalability, or performance claim. Restart is not exercised
by these refinement experiments. The rounding allowance is an experiment-specific
budget, not a universal floating-point theorem. See the
[research log](docs/research-log.md) for parameters, environment, results, and commands.

## Recorded checkpoint — Stage 5A: FDFD + GMRES integration

Committed as `69c27752fc022a04d200fa4f176dc9011db65b59`
(`feat: integrate Helmholtz FDFD with restarted GMRES`).

- Composed existing Grid2D, FDFD assembly, complex CSR matvec, manufactured
  solutions, and restarted GMRES through public APIs.
- Added four deterministic integration CTests: small, medium, altered, and
  restarted. The total suite at that checkpoint contained nine tests.
- Used algebraic forcing b=A_h*u_exact and independently checked solver status,
  true relative residual, and relative solution error on nonresonant systems.
- Exercised an actual GMRES(2) restart with a three-mode algebraic case.
- Recorded methods and results in the research log and test-coverage matrix.

This checkpoint verifies small algebraic systems; it does not establish spatial
convergence or provide a user-facing integrated solve CLI.

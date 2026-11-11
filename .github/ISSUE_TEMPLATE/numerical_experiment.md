---
name: Numerical experiment
about: Propose or document a reproducible numerical experiment
title: "[Experiment] "
---

## Research Question

[TO BE FILLED]

## Hypothesis

[TO BE FILLED]

## Mathematical Model

- Governing equation: [TO BE FILLED]
- Domain and boundary conditions: [TO BE FILLED]
- Sign/time convention: [TO BE FILLED]
- Exact or manufactured solution: [TO BE FILLED]
- Continuous vs algebraic formulation/RHS: [TO BE FILLED]

## Numerical Method

- Discretization: [TO BE FILLED]
- Grid dimensions (`nx`, `ny`): [TO BE FILLED]
- Grid spacing (`hx`, `hy`): [TO BE FILLED]
- Wavenumber and modes (`k`, `m`, `n`): [TO BE FILLED]
- Matrix structure: [TO BE FILLED]
- Solver and preconditioner: [TO BE FILLED]
- GMRES `restart`, `max_iterations`, `tolerance`: [TO BE FILLED]

## Experimental Parameters

[TO BE FILLED]

## Computational Environment

- Operating system / WSL2: [TO BE FILLED]
- CPU architecture and model: [TO BE FILLED]
- Compiler and version: [TO BE FILLED]
- CMake/build type/tool versions: [TO BE FILLED]
- Floating-point precision: [TO BE FILLED]

## Metrics

- L2 error definition: [TO BE FILLED]
- L-infinity error definition: [TO BE FILLED]
- True residual norm and normalization: [TO BE FILLED]
- Estimated residual: [TO BE FILLED]
- Final true residual: [TO BE FILLED]
- Residual history: [TO BE FILLED]
- Convergence status and iteration count: [TO BE FILLED]
- Expected convergence order: [TO BE FILLED]
- Near-resonance distance to discrete eigenvalues: [TO BE FILLED]
- Orthogonality diagnostics, if measured: [TO BE FILLED]

## Verification Strategy

[TO BE FILLED]

## Acceptance Criteria

- [ ] [TO BE FILLED]
- [ ] [TO BE FILLED]

## Expected Deliverables

[TO BE FILLED]

## Reproducibility Checklist

- [ ] Exact Git revision recorded.
- [ ] All model, grid, RHS, and solver parameters recorded.
- [ ] Computational environment recorded.
- [ ] Exact reproduction commands recorded.
- [ ] Raw data and logs preserved.
- [ ] Error and residual metrics defined.
- [ ] Acceptance criteria were fixed before inspecting results.

## Related Research Log Entry

[TO BE FILLED]

Use the reusable [Numerical Experiment Template](../../docs/templates/numerical-experiment.md)
for the full recording structure.

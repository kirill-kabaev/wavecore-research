# Numerical Experiment Template

Use this template for a reproducible WaveCore Research experiment. Replace
every `[TO BE FILLED]` marker with recorded information, or leave the marker
and explain why the field is not applicable. Do not enter values from memory or
from an unrecorded run.

## Experiment Identification

- **Experiment ID:** [TO BE FILLED]
- **Title:** [TO BE FILLED]
- **Date:** [TO BE FILLED] (use an explicit date or write `Not recorded`)
- **Researcher:** [TO BE FILLED]
- **Git commit hash:** [TO BE FILLED]
- **Related research-log entry:** [TO BE FILLED]

Record the exact source revision used for inputs, executable, and scripts. A
working-tree diff should be recorded when the revision is not clean.

## Research Question

- **Research objective:** [TO BE FILLED]
- **Hypothesis:** [TO BE FILLED]
- **Expected outcome:** [TO BE FILLED]
- **Falsification criteria:** [TO BE FILLED]

State what would count as evidence against the hypothesis before inspecting the
results. Keep expected outcomes and measured outcomes separate.

## Mathematical Model

- **Governing equation:** [TO BE FILLED]
- **Domain and boundary conditions:** [TO BE FILLED]
- **Time/sign convention:** [TO BE FILLED]
- **Exact or manufactured solution:** [TO BE FILLED]
- **Continuous formulation:** [TO BE FILLED]
- **Discrete formulation:** [TO BE FILLED]

For a 2D Helmholtz FDFD experiment, record the domain, homogeneous Dirichlet
interpretation if applicable, and the sign convention consistently. State
whether the right-hand side is continuous or algebraic. For manufactured modes,
record the mode indices `m` and `n`, wavenumber `k`, amplitude, and the
continuous/discrete eigenvalue used to form the right-hand side.

## Numerical Method

- **Discretization:** [TO BE FILLED]
- **Grid dimensions:** `nx = [TO BE FILLED]`, `ny = [TO BE FILLED]`
- **Grid spacing:** `hx = [TO BE FILLED]`, `hy = [TO BE FILLED]`
- **Physical lengths:** `Lx = [TO BE FILLED]`, `Ly = [TO BE FILLED]`
- **Wavenumber:** `k = [TO BE FILLED]`
- **Mode indices:** `m = [TO BE FILLED]`, `n = [TO BE FILLED]`
- **RHS type:** `[TO BE FILLED]` (`continuous` or `algebraic`)
- **Matrix structure:** [TO BE FILLED]
- **Solver:** [TO BE FILLED]
- **Solver parameters:** [TO BE FILLED]
- **Stopping criteria:** [TO BE FILLED]
- **Preconditioner:** [TO BE FILLED] (write `None` if none is used)

For FDFD, state whether `nx` and `ny` count interior unknowns or intervals,
and record the resulting matrix dimensions and nonzero count. For the current
five-point operator, document the x-fastest index `p = j * nx + i` when used.

For GMRES, record `restart`, `max_iterations`, `tolerance`, initial guess,
convergence status, final true residual, and the complete residual history.
Record orthogonality diagnostics such as `||V^H V - I||_F` only when they were
actually measured, including their aggregation rule and precision.

## Computational Environment

- **Operating system:** [TO BE FILLED]
- **CPU architecture and model:** [TO BE FILLED]
- **Compiler and version:** [TO BE FILLED]
- **CMake version:** [TO BE FILLED]
- **Build type:** [TO BE FILLED]
- **Compiler flags:** [TO BE FILLED]
- **Floating-point precision:** [TO BE FILLED]
- **Relevant library/tool versions:** [TO BE FILLED]

Include platform details needed to explain floating-point differences. Record
the production precision separately from any higher-precision verification
oracle.

## Experimental Procedure

- **Input generation:** [TO BE FILLED]
- **Execution steps:** [TO BE FILLED]
- **Number of repetitions:** [TO BE FILLED]
- **Random seed(s):** [TO BE FILLED] (write `Not applicable` for deterministic runs)
- **Data collection procedure:** [TO BE FILLED]
- **Reproduction commands:** [TO BE FILLED]

List commands exactly, including build options and executable paths. Describe
where raw stdout, CSV files, logs, or other artifacts are saved. Never replace
an omitted command with an assumed default.

## Verification and Validation

- **Reference solution:** [TO BE FILLED]
- **L2 error:** [TO BE FILLED]
- **L-infinity error:** [TO BE FILLED]
- **True residual norm:** [TO BE FILLED]
- **Estimated residual:** [TO BE FILLED]
- **Residual definition and normalization:** [TO BE FILLED]
- **Grid-convergence criterion:** [TO BE FILLED]
- **Expected convergence order:** [TO BE FILLED]
- **Near-resonance distance to discrete eigenvalues:** [TO BE FILLED]
- **Numerical tolerances:** [TO BE FILLED]
- **Predefined acceptance criteria:** [TO BE FILLED]

Define L2 and L-infinity norms, normalization, and whether the residual is
absolute or relative. For Helmholtz tests, record distance from `k^2` to the
relevant discrete eigenvalues. Acceptance criteria must be written before the
measured results and must not be changed after seeing them.

## Results

### Raw results table

| Case | Parameters | Status | L2 error | L-infinity error | True residual | Estimated residual | Iterations |
|---|---|---|---:|---:|---:|---:|---:|
| [TO BE FILLED] | [TO BE FILLED] | [TO BE FILLED] | [TO BE FILLED] | [TO BE FILLED] | [TO BE FILLED] | [TO BE FILLED] | [TO BE FILLED] |

- **Derived metrics:** [TO BE FILLED]
- **Plots and figure references:** [TO BE FILLED] (write `Not produced` if none)
- **Git revision associated with results:** [TO BE FILLED]

This section contains obtained results only. Do not place expected values or
acceptance criteria in the results table. Preserve raw data separately from
rounded presentation values.

## Interpretation

- **Observations:** [TO BE FILLED]
- **Interpretation:** [TO BE FILLED]
- **Alternative explanations:** [TO BE FILLED]
- **Hypothesis status:** [TO BE FILLED] (`supported`, `not supported`, or `inconclusive`)

Distinguish direct observations from explanations. A small residual is not by
itself evidence of a small solution error for an ill-conditioned system.

## Limitations

- **Numerical limitations:** [TO BE FILLED]
- **Methodological limitations:** [TO BE FILLED]
- **Environment limitations:** [TO BE FILLED]

Record missing precision, possible cancellation, conditioning, finite-size
effects, untested parameter ranges, platform dependence, and any difference
between the tested code and the intended production workflow.

## Reproducibility Checklist

- [ ] Exact code revision recorded.
- [ ] Working-tree state or relevant diff recorded.
- [ ] All numerical parameters recorded (`nx`, `ny`, `hx`, `hy`, `k`, `m`, `n`,
      RHS type, solver options, and tolerances).
- [ ] Computational environment recorded.
- [ ] Exact commands recorded.
- [ ] Raw data and logs preserved or explicitly marked unavailable.
- [ ] Error and residual metrics defined.
- [ ] Acceptance criteria were specified before results.
- [ ] Results can be independently recomputed from the recorded inputs.

## Follow-up Work

- **Open questions:** [TO BE FILLED]
- **Next experiments:** [TO BE FILLED]
- **Related issues or publications:** [TO BE FILLED]

Candidate follow-up topics include FDFD+GMRES integration, spatial refinement,
near-resonance conditioning, preconditioning, larger sparse systems, and
performance/scalability studies.

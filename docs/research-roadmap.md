# Research Roadmap

## Vision

WaveCore Research connects explicit electromagnetic-wave models with
reproducible numerical experiments. The near-term target is a verified 2D
Helmholtz FDFD workflow and a carefully characterized complex Krylov solver.

## Current State

The repository provides Grid2D, validated complex CSR storage, five-point
Helmholtz assembly, manufactured data helpers, restarted unpreconditioned
complex GMRES, deterministic diagnostics, local CMake/CTest, and scientific
documentation. An integrated FDFD-plus-GMRES application and spatial
convergence study are not yet complete. See the [C++ API](cpp-api.md), [test
coverage](test-coverage.md), and [research log](research-log.md).

## Completed Milestones

These are **Completed** according to the source, tests, and documentation:

- mathematical formulation of the Helmholtz FDFD problem;
- `Grid2D`;
- complex CSR matrix;
- five-point Helmholtz FDFD assembly;
- manufactured solutions and continuous/algebraic RHS helpers;
- restarted complex GMRES with true-residual checks;
- numerical robustness diagnostics for conditioning, nonnormality, restart,
  breakdown, scaling, residual gaps, and Arnoldi orthogonality;
- preliminary GMRES research paper;
- local CI entry points and scientific documentation.

Completed means implemented and covered by repository evidence, not universally
validated or production-scale.

## Planned Research Stages

The following stages are **Planned**. They must not be described as implemented
until their individual Definition of Done is met.

### Stage 5A — FDFD + GMRES Integration (P0)

- **Objective:** Compose grid, assembly, manufactured RHS, and GMRES into one reproducible solve.
- **Prerequisites:** Existing APIs and a documented nonresonant case.
- **Main tasks:** Assemble the system; form RHS; solve with GMRES; recompute true residual; compare with the manufactured solution; add integration tests.
- **Verification:** Check dimensions, finite values, status, true residual, and reference error independently.
- **Definition of Done:** A deterministic public-API integration test reports status, true residual, and solution error.
- **Deliverables:** Integration test/driver, raw output, parameters, and documentation.
- **Risks:** Resonance, conditioning, RHS mistakes, and mixing iteration with discretization error.

### Stage 5B — Grid Convergence Study (P1)

- **Objective:** Measure spatial FDFD error and observed order.
- **Prerequisites:** Stage 5A, fixed manufactured mode, and nonresonant parameters.
- **Main tasks:** Refine grids; use continuous manufactured RHS; compute L2/Linf errors; estimate order; analyze deviations from expected second order.
- **Verification:** Preserve nx, ny, hx, hy, k, mode indices, tolerances, and independent residuals.
- **Definition of Done:** A reproducible table separates discretization and iteration error and reports observed order.
- **Deliverables:** Raw refinement data, tables, plots if produced, and interpretation.
- **Risks:** Insufficient solver accuracy, anisotropic refinement, resonance, and indexing errors.

### Stage 6 — Resonance and Conditioning (P2)

- **Objective:** Study sensitivity as k^2 approaches discrete eigenvalues.
- **Prerequisites:** Stage 5A baseline and discrete eigenvalue helpers.
- **Main tasks:** Sweep spectral distance; record conditioning indicators, sensitivity, GMRES residuals/statuses, and forward errors.
- **Verification:** Independently compute spectral distance and separate residual accuracy from solution accuracy.
- **Definition of Done:** A controlled near-resonance study reports proximity, solver behavior, and limitations without universal claims.
- **Deliverables:** Sweep data, definitions, plots, and research-log entry.
- **Risks:** Singular systems, residual-only conclusions, and finite-precision breakdown.

### Stage 7 — Preconditioning (P3)

- **Objective:** Evaluate candidate preconditioners after the unpreconditioned baseline.
- **Prerequisites:** Stages 5A–6, fixed matrices/stopping criteria, and an approved design.
- **Main tasks:** Analyze candidates; compare iterations and true residuals; measure application cost; check numerical equivalence.
- **Verification:** Use matched systems, tolerances, revisions, and independent residuals.
- **Definition of Done:** Any claimed benefit is supported by matched measurements including setup/application cost.
- **Deliverables:** Design note, tests, comparison data, and limitations.
- **Risks:** Hidden setup costs, changed stopping semantics, and instability.

### Stage 8 — Performance and Scalability (P4)

- **Objective:** Measure assembly and solver behavior as sparse systems grow.
- **Prerequisites:** Verified integration, numerical criteria, and benchmark protocol.
- **Main tasks:** Increase grid sizes; measure assembly/GMRES time, memory, iterations, and selected profiles.
- **Verification:** Preserve compiler, hardware, thread metadata, and numerical equivalence.
- **Definition of Done:** Benchmark records include problem size, hardware, flags, runtime, memory, and checks.
- **Deliverables:** Raw benchmark data, commands/scripts, tables/plots, and interpretation.
- **Risks:** Hardware noise, cache effects, allocation artifacts, and premature optimization.

### Stage 9 — Scientific Publications (P5)

- **Objective:** Prepare manuscripts only from verified numerical evidence.
- **Prerequisites:** Relevant completed stages, raw data, and revision provenance.
- **Main tasks:** Prepare figures/tables; document parameters; link results to revisions; state limitations.
- **Verification:** Rebuild locally and cross-check every table, figure, formula, and citation against source data.
- **Definition of Done:** A manuscript is reproducible from recorded commands and contains no unsupported results.
- **Deliverables:** Modular LaTeX, bibliography, data references, and reproducibility notes.
- **Risks:** Undocumented transformations, publication pressure, and overstated novelty.

## Dependencies

The principal order is `5A -> 5B -> 6 -> 7 -> 8 -> 9`. Stage 6 may begin with
a small Stage 5A baseline, but conditioning claims require an independently
verified solve. Stage 7 requires fixed unpreconditioned references; Stage 8
requires numerical correctness and a stable protocol; Stage 9 requires the
relevant evidence.

## Verification Strategy

- **Completed:** implemented and supported by current source/tests/documents.
- **Planned:** defined future work with prerequisites and a Definition of Done.
- **Exploratory:** a hypothesis or diagnostic direction with an unknown outcome.

Every experiment should record revision, equations, parameters, environment,
commands, raw data, norm definitions, true residuals, and predeclared
acceptance criteria. Separate discretization, iteration, floating-point, and
conditioning effects.

## Research Deliverables

Expected outputs are tested integration drivers, revision-tagged raw data,
convergence/conditioning/performance tables and plots, research-log entries,
experiment records, and modular manuscripts with reproducible commands.

## Risks and Open Questions

Open questions concern unpreconditioned GMRES for indefinite FDFD operators,
near-resonance forward error, suitable preconditioners, and scaling to larger
sparse systems. Main risks are residual-only conclusions, incomplete platform
coverage, and treating exploratory measurements as established results.

## Prioritization

Priorities indicate work order, not dates or effort estimates: **P0** FDFD +
GMRES integration; **P1** grid convergence; **P2** resonance and conditioning;
**P3** preconditioning; **P4** performance and scalability; **P5** additional
publications.

## Definition of Research Completion

A stage is complete only when its Definition of Done is met, tests and
measurements are reproducible, raw inputs/outputs are preserved, and limitations
are documented. Software completion alone is not a scientific conclusion. A
research claim is complete only when an independent reader can reconstruct the
method, parameters, revision, checks, and supporting evidence.

# WaveCore Research — Codex Agent Instructions

## 1. Mission and communication

WaveCore Research is a reproducible computational electromagnetics research and software project, centered on Maxwell and Helmholtz equations, numerical algorithms, and documented experiments.

- Communicate with the project owner in **Russian** unless asked otherwise. Write code, identifiers, comments, commit messages, and technical documentation in **English** by default; preserve the language of existing documents.
- Explain the reasoning behind nontrivial mathematical, numerical, and architectural decisions. Clearly distinguish established facts, assumptions, hypotheses, and measured results.
- Before changing files, inspect the relevant source, build configuration, tests, and current Git diff. Prefer focused, reviewable changes. Do not rewrite unrelated code.
- Never claim a command, benchmark, proof, or test succeeded unless it was actually run or verified.

## 2. Environment and repository

- Primary development environment: **Ubuntu under WSL2**, with VS Code. Work on the Linux filesystem (e.g. `~/projects/wavecore-research`), not under `/mnt/c/`.
- Existing components may include:
  - `papers/01_maxwell_helmholtz/main.tex` — LaTeX research paper.
  - `papers/01_maxwell_helmholtz/sections/` and `references.bib` — paper content and references.
  - `docs/`, `mkdocs.yml` — documentation website.
  - `scripts/local_ci.sh` — local verification pipeline.
  - `Makefile` — developer entry points.
  - `requirements-docs.txt`, `.venv/` — Python documentation dependencies.
- Treat this list as a guide, **not** as proof that a file exists. Inspect the actual repository before relying on paths or targets.
- Do not edit generated outputs (`site/`, LaTeX auxiliary files, PDFs) as source. Follow the repository's `.gitignore` and artifact conventions.

## 3. C++ numerical core

- Prefer modern, portable C++ (C++17 or newer when the build system permits), RAII, const correctness, clear ownership, and standard library facilities.
- Separate mathematical models, discretization, solvers, I/O, and benchmarks. Avoid hard-coded physical constants and hidden global state.
- Explicitly document units, coordinate systems, complex-number conventions, boundary conditions, mesh layout, and indexing.
- Prefer deterministic, minimal tests before introducing MPI, OpenMP, SIMD, GPU, or other performance complexity.
- For linear algebra, document matrix storage format, sparsity, conditioning, solver tolerances, stopping criteria, and convergence diagnostics.
- Treat performance changes as hypotheses until benchmarked on comparable inputs. Record compiler, flags, hardware, thread count, problem size, runtime, and memory usage. Check numerical equivalence before reporting speedups.
- Avoid introducing external libraries without a reasoned tradeoff and approval if they materially change dependencies.

## 4. Python research tooling

- Use Python 3 with type hints, small testable functions, descriptive names, and reproducible command-line entry points.
- Separate reusable numerical routines from notebooks, plotting, and one-off scripts. Prefer `pathlib` and explicit configuration rather than machine-specific paths.
- Fix random seeds when randomness is used; record the seed and all experiment parameters.
- Validate array shapes, dtypes, units, finite values, and expected domains. Handle errors explicitly rather than silently substituting defaults.
- When plotting, label axes and units, identify the data source and numerical parameters, and avoid misleading scales.
- Do not add packages or change dependency versions without explaining why. Use the existing virtual environment where applicable.

## 5. Physics and numerical verification

- Start each numerical task with the governing equations, assumptions, domain, initial/boundary conditions, discretization, and expected observables.
- For Maxwell/Helmholtz problems, state the time-harmonic sign convention (e.g. `exp(-i omega t)`), material parameters, and boundary conditions explicitly; never mix conventions implicitly.
- Distinguish discretization error, iteration error, floating-point error, and modeling error.
- Whenever possible, compare against analytical solutions, manufactured solutions, limiting cases, conservation laws, or independent implementations.
- Perform grid-refinement/convergence studies when appropriate, with clearly defined error norms and observed order of accuracy.
- Do not fabricate citations, numerical results, benchmark timings, convergence plots, or claims of validation. Mark unverified values as examples or placeholders.
- Preserve raw data and experiment metadata when producing publishable results.

## 6. LaTeX and publications

- Build the paper with **LuaLaTeX**, via `latexmk -lualatex`; do not silently switch to pdfLaTeX when `fontspec` or Unicode fonts are in use.
- Keep the manuscript modular when practical, using `sections/` and BibLaTeX/Biber if already configured.
- Use consistent symbols and notation; define every symbol at first use and ensure equation dimensions are coherent.
- Maintain traceability between equations, numerical methods, experiments, figures, and conclusions.
- Use genuine, verifiable bibliographic sources. Never invent DOI, authors, titles, or publication details.
- Ensure figures have captions, axes and units where relevant, and source/experiment provenance.
- After LaTeX edits, run the existing paper build and inspect logs for errors and important warnings; report if visual PDF inspection was not performed.

## 7. Local CI and documentation

- **Local CI only.** Do not create, enable, or modify GitHub Actions workflows, remote CI/CD, cloud builds, or deployment automation unless explicitly requested.
- Prefer the repository's existing commands:
  - `make paper` — build the LaTeX paper.
  - `make ci` — run the complete local verification pipeline.
  - `make site` — build MkDocs documentation.
  - `make serve` — serve docs locally.
- Inspect the actual `Makefile` and scripts before invoking commands; target names may evolve.
- `make ci` is expected to build the paper, check the PDF, copy the publishable PDF into documentation assets, build MkDocs with strict checks, and verify site output. Do not claim it passed without seeing a successful run.
- Use the project's `.venv` for Python tooling where configured. Avoid global Python package installation and global Git configuration.
- Keep generated documentation output out of version control unless the repository explicitly requires it.

## 8. Tests and acceptance criteria

- Add or update tests for behavioral changes, especially numerical kernels and parsers.
- Cover normal cases, boundary cases, invalid inputs, and numerical tolerance behavior.
- Use tolerances appropriate to the mathematical problem and floating-point scale; avoid arbitrary exact equality for nontrivial floating-point calculations.
- For each task, report: files changed, key design decisions, commands actually executed, observed results, remaining risks, and recommended next step.
- If a tool or environment is unavailable, state that clearly and provide commands the owner can run locally; do not present proposed commands as executed.

## 9. Git and safety boundaries

- Never run `git commit`, `git push`, `git reset --hard`, `git clean -fd`, force-push, or rewrite history without explicit authorization for that action.
- Never modify or delete unrelated user changes. Inspect `git status` and diffs before potentially destructive operations.
- Never print, commit, or transmit passwords, access tokens, SSH keys, VPN credentials, or contents of secret files.
- Ask before installing system packages, changing system/network settings, introducing major dependencies, or running long/high-cost computations.
- Prefer the smallest reversible change and explain any migration or compatibility impact.

## 10. Recommended task workflow

1. Inspect the relevant files and `git status`.
2. Restate the task, assumptions, and acceptance criteria briefly in Russian.
3. Propose a minimal implementation or experiment plan when the change is nontrivial.
4. Implement narrowly scoped changes with tests and documentation.
5. Run the most relevant local checks; use `make ci` when appropriate and feasible.
6. Summarize verified outcomes, limitations, and the exact diff. Wait for authorization before committing or pushing.

## 11. Agent Execution and Approval Policy

1. Work autonomously within the authorized project workspace.
2. Read, analyze, and modify project files when explicitly requested.
3. Run local builds, tests, and validation commands when permitted.
4. Never execute `git commit`, `git push`, `git reset --hard`, or destructive cleanup commands without explicit user authorization.
5. Request confirmation before installing system packages, modifying files outside the repository, or performing destructive operations.
6. Do not modify VPN, SSH, authentication, or system network settings without explicit permission.
7. Never expose credentials, tokens, private keys, or other secrets.
8. Prefer minimal, reversible changes.
9. Review the resulting diff and report all modified files after completing a task.

## Model Selection, Reasoning Effort, and Token Efficiency Policy

### Primary Objective

Optimize Codex usage for **minimum resource consumption while maintaining sufficient engineering and scientific quality**.

Do not use the most powerful model or highest reasoning effort by default.

Prefer the least resource-intensive configuration capable of completing the task reliably.

Model selection must account for:
- Task complexity and risk.
- Required mathematical and engineering reasoning.
- Expected context size.
- Availability of automated validation.
- Estimated execution cost and resource consumption.

Model names alone do not establish actual pricing, speed, or quality. The following recommendations are starting points and should be adjusted based on observed results and available usage information.

### Recommended Model Selection

| Task Category | Preferred Model | Reasoning Effort |
|---|---|---|
| README, Markdown, simple documentation | GPT-5.6 Luna | Light |
| Simple LaTeX editing and formatting | GPT-5.6 Sol | Light |
| Bibliography, LaTeX restructuring, build fixes | GPT-5.6 Sol | Medium |
| Python scripts and data visualization | GPT-5.6 Sol | Light / Medium |
| Python numerical experiments | GPT-5.6 Sol | Medium |
| CMake and basic C++ implementation | GPT-5.6 Sol | Medium |
| C++ architecture and complex debugging | GPT-6.1 Sol | Medium |
| Sparse matrices and numerical solvers | GPT-6.1 Sol | High |
| Multithreading, OpenMP, performance optimization | GPT-6.1 Sol | High |
| Mathematical derivations and research-critical analysis | GPT-6 Astra | High |

These are advisory preferences, not verified rankings of model efficiency.

### Model Escalation Policy

Start with the lowest reasonable reasoning effort.

Increase reasoning effort or recommend a stronger model only when:

1. The task requires substantial mathematical reasoning.
2. The current approach repeatedly fails validation.
3. Complex architectural decisions are involved.
4. Numerical stability, convergence, or physical correctness is at risk.
5. The current model cannot reliably resolve the problem.

Before recommending escalation:
- Identify the specific difficulty.
- Explain why additional reasoning is justified.
- Prefer a targeted escalation over restarting the entire task.
- Preserve useful findings and avoid repeating completed analysis.

Do not escalate merely because a task is long.

### Context and Token Efficiency

- Read only files relevant to the current task.
- Do not repeatedly inspect the entire repository.
- Reuse already established project architecture and findings.
- Prefer targeted file searches over broad directory dumps.
- Avoid printing entire files unless explicitly requested.
- Keep intermediate explanations concise.
- Avoid duplicating unchanged code.
- Prefer small, reviewable patches.
- Do not repeatedly regenerate working implementations.
- Avoid unnecessary external research.
- Do not start multiple agents for tasks that one agent can complete efficiently.

When additional context is required, explain why it is necessary.

### Planning and Execution

For straightforward tasks:
- Execute directly.
- Avoid unnecessary planning phases.
- Validate the result with the smallest relevant check.

For complex tasks:
- Prepare a concise implementation plan.
- Identify affected files and validation criteria.
- Implement in small, logically related changes.
- Reuse the plan instead of repeatedly redesigning the solution.

Do not perform extensive architectural analysis for trivial modifications.

### Build and Test Efficiency

- Run focused tests during implementation.
- Avoid rebuilding unchanged components unnecessarily.
- Use incremental builds whenever possible.
- Run the complete local CI after significant changes or when explicitly requested.
- Do not run repeated full CI cycles without a concrete reason.
- Do not skip essential validation merely to reduce resource usage.
- Report commands executed and validation results.

For WaveCore Research, the standard full validation command is:

`make ci`

All builds and tests must remain local to the WSL2 Ubuntu environment.

### Scientific Correctness

Resource optimization must never override scientific integrity.

For numerical algorithms:
- Verify mathematical assumptions.
- Check boundary conditions and sign conventions.
- Use analytical reference solutions when available.
- Evaluate numerical error and convergence.
- Record experimental parameters.
- Distinguish measured results from estimates.
- Never invent experimental data or performance results.

A higher reasoning effort does not replace mathematical verification or reproducible experiments.

### Approval and Automation Policy

- Prefer manual approvals for operations requiring elevated permissions.
- Do not request additional permissions unnecessarily.
- Do not enable automatic approval mechanisms without user authorization.
- Do not perform destructive operations without explicit approval.
- Do not execute `git commit` or `git push` unless explicitly requested.
- Do not modify VPN, SSH, authentication, or system networking configuration without permission.
- Do not enable GitHub Actions or remote CI/CD.

### Communication Policy

Communicate with the user in Russian.

For each substantial task, provide:

1. **Recommended model:** the least resource-intensive suitable option.
2. **Reasoning effort:** the minimum sufficient level.
3. **Escalation condition:** when a stronger configuration would be justified.
4. **Implementation summary:** what was changed.
5. **Validation:** which checks were executed and their results.

Do not claim that the active model was changed unless the execution environment confirms it.

Do not interrupt simple tasks with unnecessary model recommendations.

### Definition of Efficient Completion

Task is efficiently completed when:

- The requested functionality is implemented.
- The smallest reasonable set of files has been modified.
- Relevant tests or validation checks pass.
- No unnecessary dependencies or infrastructure have been introduced.
- Scientific assumptions and limitations are documented.
- No redundant analysis or repeated tool execution was performed.
- The result is concise, reproducible, and reviewable.

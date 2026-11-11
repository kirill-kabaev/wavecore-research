# WaveCore Research

## Overview

WaveCore Research is a reproducible computational electromagnetics project. The
current numerical core focuses on two-dimensional Helmholtz discretization and
iterative solution of complex-valued sparse linear systems. Source code,
verification tests, and research manuscripts are maintained in one repository.

## Research Objectives

The project develops a transparent path from mathematical specifications to
tested numerical building blocks:

- model Maxwell/Helmholtz wave problems with explicit conventions and units;
- implement small, portable C++ numerical kernels;
- verify discretization and solver behavior with deterministic tests;
- record methods and limitations in reproducible research documents.

End-to-end FDFD-plus-GMRES integration and a systematic spatial-convergence
study are planned work. They are not presented as current functionality.

## Current Implementation

The C++20 library currently provides:

- `Grid2D`, a validated uniform rectangular grid containing interior unknowns,
  with x-fastest indexing `p = j * nx + i`;
- an immutable complex-valued CSR matrix with strict structural validation and
  sparse matrix-vector multiplication;
- assembly of the five-point finite-difference Helmholtz operator
  `A = -Delta_h - k^2 I` for homogeneous Dirichlet boundaries;
- continuous and discrete modal eigenvalue helpers and manufactured right-hand
  sides for the Helmholtz operator;
- sequential, unpreconditioned restarted complex GMRES with Hermitian
  Arnoldi orthogonalization, two-pass modified Gram--Schmidt, complex Givens
  rotations, true-residual checks, and diagnostic status information;
- deterministic numerical verification tests covering the grid, CSR matrix,
  operator assembly, and GMRES behavior;
- a local CMake/CTest build and test infrastructure.

The command-line application is currently a minimal project-information entry
point; it does not yet run an integrated Helmholtz solve.

## Repository Structure

- [`CMakeLists.txt`](CMakeLists.txt) — C++20 targets, warnings, and CTest registration.
- [`include/wavecore/helmholtz/`](include/wavecore/helmholtz/) — public numerical APIs.
- `src/helmholtz/` — library implementations.
- `src/apps/` — command-line applications.
- `tests/helmholtz/` — deterministic C++ verification tests.
- [`experiments/helmholtz_2d/README.md`](experiments/helmholtz_2d/README.md) — mathematical specification of the first 2D experiment.
- `papers/01_maxwell_helmholtz/` — foundational Maxwell/Helmholtz manuscript.
- `papers/02_gmres_numerical_observations/` — preliminary GMRES numerical-observations manuscript.
- `docs/` and [`mkdocs.yml`](mkdocs.yml) — the local documentation website.
- [`Makefile`](Makefile) — local paper, CI, and documentation entry points.

## Requirements

- CMake 3.20 or newer and a C++20 compiler (GCC or Clang are supported by the
  warning configuration).
- Make and a standard Unix shell.
- LuaLaTeX, `latexmk`, and Biber for the manuscripts.
- Python 3 with MkDocs and the configured documentation plugins in the
  repository's `.venv` environment for website builds.

No third-party C++ library is required by the current implementation.

## Building the C++ Project

Use an out-of-source build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel 2
```

The build produces the `wavecore` library and the `helmholtz_2d` executable.

## Running Tests

After configuring and building, run:

```bash
ctest --test-dir build --output-on-failure
```

The registered tests are `wavecore.library_metadata`, `wavecore.grid`,
`wavecore.csr_matrix`, `wavecore.assembly`, and `wavecore.gmres`.

## Building Scientific Papers

The first manuscript is built through the Makefile:

```bash
make paper
```

To build the second manuscript directly with LuaLaTeX and its BibLaTeX/Biber
references:

```bash
cd papers/02_gmres_numerical_observations
latexmk -lualatex -interaction=nonstopmode -file-line-error -halt-on-error main.tex
cd ../..
```

Generated PDFs and auxiliary LaTeX files are build products and are excluded
from version control by the repository conventions.

## Local CI

WaveCore uses local verification only. Run the repository's complete local
pipeline with:

```bash
make ci
```

The pipeline is defined by the existing local scripts and does not use GitHub
Actions, hosted CI, or deployment automation.

## Documentation Website

Build the MkDocs site locally with the existing Makefile target:

```bash
make site
```

To serve it during development:

```bash
make serve
```

The generated `site/` directory is a local artifact and is not a source file.

## Publications

These are repository manuscripts and technical research documents; the project
does not claim that they have been published in a peer-reviewed journal.

- [Physical foundations of electromagnetic-wave numerical modeling](papers/01_maxwell_helmholtz/main.tex)
- [Numerical Robustness of Restarted GMRES for Complex-Valued Linear Systems: Preliminary Computational Experiments](papers/02_gmres_numerical_observations/main.tex)

## Project Status

The grid, complex CSR storage, five-point Helmholtz assembly, restarted GMRES,
and their deterministic verification tests are implemented. The integrated
FDFD solve driven by GMRES, larger-scale performance studies, and spatial
convergence experiments remain planned extensions.

## Reproducibility

Start from a clean checkout, use the commands above, and retain the compiler,
CMake version, build type, platform, and numerical parameters with any new
measurement. The tests are deterministic and do not require network access.
Separate discretization error, iteration error, floating-point effects, and
modeling assumptions when interpreting results. Do not treat a passing unit
test as evidence of universal solver stability or of a completed scientific
validation study.

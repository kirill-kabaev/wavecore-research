# Developer Commands

This page is the authoritative quick reference for commands exposed by the
current repository configuration. Run commands from the repository root unless
a different working directory is shown. Commands below are documented from
`Makefile`, `CMakeLists.txt`, `mkdocs.yml`, and the existing paper build setup;
they are not claims that the commands have been executed in this checkpoint.

## Environment and Prerequisites

WaveCore targets Ubuntu/WSL2 development on the Linux filesystem. The current
configuration requires:

- CMake 3.20 or newer;
- a C++20 compiler (GCC or Clang receive `-Wall -Wextra -Wpedantic`);
- GNU Make and a POSIX shell;
- LuaLaTeX, `latexmk`, and Biber for BibLaTeX manuscripts;
- Python 3 with MkDocs Material and configured plugins in `.venv`.

From WSL2, keep the checkout under a Linux path such as
`~/projects/wavecore-research`, open it in VS Code with the WSL integration,
and use a WSL terminal. This avoids the slower and sometimes different file
semantics of `/mnt/c`. Build products belong in the repository's out-of-source
`build/` directory; MkDocs writes its generated site to `site/`.

## CMake Configuration

**Working directory:** repository root.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
```

This configures an out-of-source CMake build. The project requires CMake 3.20,
uses C++20, rejects in-source builds, and creates the `wavecore` static library,
the `helmholtz_2d` executable, and CTest targets when `BUILD_TESTING` is on.
Configuration creates or updates files under `build/`; it does not modify
source files. A common error is running CMake from the wrong directory or
attempting `cmake .`, which is explicitly unsupported.

## Building the C++ Project

**Working directory:** repository root, after configuration.

```bash
cmake --build build --parallel 2
```

This compiles the library, application, and registered test executables. It
creates build artifacts under `build/`. The command does not run tests. Typical
failures are a missing C++20 compiler, a stale or incorrectly configured build
directory, or compiler diagnostics treated as errors by a local toolchain.

## Running CTest

**Working directory:** repository root; the configured build directory is
`build/`.

```bash
ctest --test-dir build --output-on-failure
```

The current CMake registration contains five tests: library metadata, Grid2D,
CSR matrix, Helmholtz assembly, and GMRES. CTest executes the already-built
test binaries and reports failures; it does not rebuild the project. Configure
and build first if `build/` is absent or stale.

## Building Scientific Papers

### First manuscript

**Working directory:** repository root for the Make target.

```bash
make paper
```

The Makefile changes into `papers/01_maxwell_helmholtz/` and runs:

```bash
latexmk -lualatex -interaction=nonstopmode \
  -file-line-error -halt-on-error main.tex
```

LuaLaTeX compiles the modular manuscript. Because the paper uses BibLaTeX with
the Biber backend, `latexmk` invokes Biber as part of the normal dependency
cycle when bibliography data require it; there is no separate Make target for
Biber. PDFs and auxiliary LaTeX files are generated artifacts.

### Second manuscript

**Working directory:** `papers/02_gmres_numerical_observations/`.

```bash
latexmk -lualatex -interaction=nonstopmode \
  -file-line-error -halt-on-error main.tex
```

This is the direct LuaLaTeX build for the GMRES manuscript. Its BibLaTeX/Biber
backend is handled by `latexmk`. The repository does not provide a dedicated
root Make target for this second paper.

Missing LuaLaTeX, `latexmk`, Biber, fonts, or bibliography data are typical
failure causes. Paper builds can take longer than metadata or documentation
checks and write PDFs and auxiliary files beside the manuscript sources.

## Building MkDocs

**Working directory:** repository root.

The explicit Make target is:

```bash
make site
```

It runs `.venv/bin/python -m mkdocs build --strict`, producing the generated
website in `site/`. The equivalent direct command is:

```bash
.venv/bin/python -m mkdocs build --strict
```

The strict mode treats broken documentation references and configuration issues
as failures. The repository configuration explicitly includes the numerical
experiment template under `docs/templates/`, because MkDocs otherwise excludes
that directory by default. This command does not build C++ or papers.

For an interactive local server, the supported Make target is:

```bash
make serve
```

It runs `.venv/bin/python -m mkdocs serve` and is intended for development, not
for producing a committed site artifact.

## Local Continuous Integration

**Working directory:** repository root.

```bash
make ci
```

This is WaveCore's local CI entry point. It prepends `.venv/bin` to `PATH` and
invokes `scripts/local_ci.sh`. The pipeline is local-only; GitHub Actions,
hosted CI, and deployment automation are not used. `make ci` can perform
multiple document/site checks and may take longer than a single documentation
build. Consult the script output for the exact checks performed by the current
revision; do not assume it replaces a separate CMake/CTest run unless the
output explicitly shows those steps.

## Cleaning Generated Files

The only supported cleaning target is:

```bash
make clean
```

It changes into `papers/01_maxwell_helmholtz/` and runs `latexmk -C`, cleaning
that paper's LaTeX auxiliary files. It does not clean the CMake `build/`
directory, the second paper's artifacts, or the MkDocs `site/` directory.

No broader repository cleanup command is defined here. Do not substitute an
unreviewed recursive deletion command; inspect generated paths and project
conventions before removing artifacts manually.

## Common Troubleshooting

- **CMake reports an in-source build:** use `cmake -S . -B build ...` from the
  repository root.
- **Build directory is missing:** configure first, then run `cmake --build`.
- **CTest cannot find tests:** rebuild the configured tree and verify that
  `BUILD_TESTING` was enabled.
- **Paper bibliography is unresolved:** confirm LuaLaTeX, `latexmk`, Biber, and
  the paper's `references.bib` are available; rerun the same `latexmk` command.
- **MkDocs cannot start:** use the repository's `.venv/bin/python` rather than
  assuming a globally installed MkDocs.
- **Strict documentation failure:** inspect the reported page, navigation
  target, or link; do not suppress strict mode to hide a broken reference.
- **WSL2 path or permission differences:** work under the Linux filesystem and
  run the command in the WSL terminal associated with the checkout.

## Command Reference Table

| Command | Working directory | Purpose / result | Artifacts and typical constraints |
|---|---|---|---|
| `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug` | Repository root | Configure out-of-source C++20 build and CTest targets | Writes `build/`; requires CMake and compiler; source-only check, no compilation |
| `cmake --build build --parallel 2` | Repository root | Build library, application, and test executables | Writes `build/`; requires configured tree and compiler; can take longer than inspection commands |
| `ctest --test-dir build --output-on-failure` | Repository root | Run registered C++ tests and print failures | Reads/runs `build/`; requires successful build; does not rebuild |
| `make paper` | Repository root | Build first paper with LuaLaTeX via `latexmk` | Creates PDF and LaTeX auxiliaries; requires LuaLaTeX, Biber, fonts |
| `latexmk -lualatex -interaction=nonstopmode -file-line-error -halt-on-error main.tex` | `papers/02_gmres_numerical_observations/` | Build second paper directly | Creates PDF and auxiliaries; requires paper toolchain and bibliography data |
| `make site` | Repository root | Run MkDocs strict build | Writes `site/`; requires `.venv` MkDocs environment |
| `make serve` | Repository root | Serve MkDocs locally for interactive review | Runs a local server; not a publication artifact |
| `make ci` | Repository root | Run the configured local CI pipeline | May build/check multiple artifacts; local-only, no GitHub Actions |
| `make clean` | Repository root | Run `latexmk -C` for the first paper | Cleans only first-paper LaTeX auxiliaries; not a general repository cleanup |

Commands that create artifacts are the CMake configure/build, paper builds,
MkDocs build, and the local CI steps they invoke. `ctest` primarily verifies
existing build outputs; `make serve` runs a development server. Always record
the command, working directory, revision, and relevant tool versions with
reproducible research results.

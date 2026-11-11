# Numerical Verification Protocol

## 1. Scope and Mathematical Model

This protocol consolidates evidence at baseline
58fcf252f97e69b60b6bf1efe2fd0ece9286c8a6. It distinguishes:

- **Mathematical deduction:** properties of the ideal uniform-grid Dirichlet stencil in exact arithmetic.
- **Numerical observation:** recorded values for specific experiments and an identified environment.
- **Automated verification:** assertions exercised by registered CTest cases.
- **Not verified:** general stability, arbitrary PDE inputs, near-resonant Helmholtz behavior, and cross-platform reproducibility.

Stage 5A verifies algebraic integration; Stage 5B verifies square-grid spatial
convergence; Stage 5C.1 tests acceptance semantics; Stage 5C.2 verifies
rectangular spatial convergence. Passing these tests is not a formal proof
of program correctness for arbitrary problems.

The scalar model is the z-invariant TM_z electric field \(u=E_z\), with
time dependence \(\exp(-i\omega t)\), in a homogeneous, isotropic, lossless medium:

$$
-\Delta u-k^2u=f\ \text{in }\Omega=(0,L_x)\times(0,L_y),\qquad
u|_{\partial\Omega}=0,\qquad k=\omega\sqrt{\mu\epsilon}>0.
$$

Here \(\omega\) is angular frequency, \(\mu\) permeability and \(\epsilon\)
permittivity; both material parameters are real positive constants.
PEC walls force the tangential electric component \(E_z\) to vanish.
Coordinates/lengths are in metres, \(k\) in \(\mathrm{m}^{-1}\),
eigenvalues in \(\mathrm{m}^{-2}\). With amplitude in V/m, \(f\) has units
V/m cubed; weighted two-dimensional L2 errors have units V and Linf errors
V/m. Residual ratios are dimensionless.

This is a forced scalar model, not a verification of a complete Maxwell
source model, heterogeneous media, radiation conditions, or PML.
Continuous Dirichlet eigenvalues are
\(\lambda_{mn}=(m\pi/L_x)^2+(n\pi/L_y)^2\), \(m,n\ge1\).
Uniqueness requires \(k^2\) outside this spectrum; proximity can amplify
errors. Discrete resonance must be checked separately.

## 2. Discrete FDFD Operator

Grid dimensions \(n_x,n_y\) count **interior nodes**, not intervals:

$$
h_x=\frac{L_x}{n_x+1},\quad h_y=\frac{L_y}{n_y+1},\quad
x_i=(i+1)h_x,\quad y_j=(j+1)h_y,\quad p=j n_x+i.
$$

Here \(0\le i<n_x\), \(0\le j<n_y\); x varies fastest.
For \(A_h=-\Delta_h-k^2I\),

$$
(A_hu)_{ij}=
\left(\frac2{h_x^2}+\frac2{h_y^2}-k^2\right)u_{ij}
-\frac{u_{i-1,j}+u_{i+1,j}}{h_x^2}
-\frac{u_{i,j-1}+u_{i,j+1}}{h_y^2}.
$$

Boundary values are zero; absent neighbors are omitted, never wrapped to the
next grid row. There are \(n_xn_y\) unknowns and

$$
\operatorname{nnz}=n_xn_y+2(n_x-1)n_y+2n_x(n_y-1).
$$

The matrix uses complex-double CSR values and size_t indices, with strictly
increasing column indices per row. The current real-coefficient operator is
real symmetric, although storage and GMRES support complex data.
CSR validation rejects malformed structure and nonfinite stored values.
See the [API](cpp-api.md) for signatures and exceptions.

For resolved modes \(1\le m\le n_x,\ 1\le n\le n_y\), discrete sine analysis gives

$$
\lambda_{mn,h}=
\frac4{h_x^2}\sin^2\frac{m\pi h_x}{2L_x}
+\frac4{h_y^2}\sin^2\frac{n\pi h_y}{2L_y},\qquad
\mu_{mn,h}=\lambda_{mn,h}-k^2.
$$

These exact-arithmetic stencil eigenvalues differ from continuous eigenvalues
and from a floating-point eigendecomposition of stored CSR coefficients.
Assembly tests compare small matrices, symmetry, counts, and modal identities.
Convergence drivers compute reference eigenvalues independently of production
eigenvalue helpers.

## 3. GMRES Verification

The solver is sequential, unpreconditioned restarted complex GMRES(m).
Arnoldi uses two passes of modified Gram–Schmidt and the Hermitian product

$$
\langle v,w\rangle=\sum_p\overline{v_p}w_p,\qquad
AV_q=V_{q+1}\overline H_q.
$$

The columns of \(V_q\) are the q Arnoldi basis vectors; \(\overline H_q\)
is the (q+1)-by-q upper Hessenberg matrix. The bar denotes the extended
Hessenberg matrix here, not entrywise conjugation.
The small least-squares problem uses complex Givens rotations. The convention is
\(G=\left[\begin{smallmatrix}c&s\\-\overline s&c\end{smallmatrix}\right]\),
with real \(c\), \(c^2+|s|^2=1\). For \(a=0,b\ne0\),
\(c=0,\ s=\overline b/|b|\), so \(G(0,b)^T=(|b|,0)^T\).
Tests exercise complex systems and the zero-first-entry case.

The stopping measure is

$$
\rho=\frac{\|b-A_h\widehat u_h\|_2}{\|b\|_2}\quad(b\ne0).
$$

For exactly zero RHS it is the absolute residual norm. The least-squares
estimate alone cannot establish convergence: the solver recomputes the true
residual. Drivers independently recompute it with a fresh CSR matvec and
wide-precision norm accumulation. The rectangle also uses a separate
complex-long-double product over stored CSR entries and checks both residuals.
This independence concerns solver bookkeeping, not an entirely separate PDE
implementation.

CONVERGED is mandatory for accepted PDE experiments; NOT_CONVERGED and
BREAKDOWN fail acceptance. Happy breakdown succeeds only after true-residual
confirmation. max_iterations is global across restart cycles; restart bounds
the Krylov dimension per cycle. See the [GMRES API](cpp-api.md).
The GMRES unit test compares a local one-pass Arnoldi control with the
production two-pass basis using \(\|V^HV-I\|_F\); this is a diagnostic,
not a universal orthogonality bound.

## 4. Manufactured Solutions

Let \(\phi_{mn}=\sin(m\pi x/L_x)\sin(n\pi y/L_y)\), with amplitude \(U_0\).

- **Algebraic forcing:** \(b=A_h(U_0\phi_{mn}|_h)\), cross-checked against
  \((\lambda_{mn,h}-k^2)U_0\phi_{mn}|_h\). The target is the sampled vector.
  This tests the linear solve, not spatial truncation error.
- **Continuous forcing:** \(f=(\lambda_{mn}-k^2)U_0\phi_{mn}\), sampled
  without using \(A_hu^*|_h\). The ideal discrete solution is

$$
u_h=\frac{\lambda_{mn}-k^2}{\lambda_{mn,h}-k^2}U_0\phi_{mn}|_h,\qquad
d_{mn,h}=\frac{\lambda_{mn}-\lambda_{mn,h}}{\lambda_{mn,h}-k^2}.
$$

Stage 5B uses the unit square, k=1 inverse metre, and either mode (1,1) or
\(u^*=\phi_{11}-0.3\phi_{23}\). With lengths numerically expressed in metres,
the mixture RHS is
\((2\pi^2-k^2)\phi_{11}-0.3(13\pi^2-k^2)\phi_{23}\);
the implicit unit-amplitude factor carries the field units.
Each mode has its own discrete scaling. Stage 5C.2 uses mode (1,1) on
Lx=1 m, Ly=1.5 m, with \(\lambda_{11}=13\pi^2/9\) in inverse square metres.

The decomposition is

$$
\widehat u_h-u^*|_h =
\underbrace{\widehat u_h-u_h}_{e_{\rm alg}}+
\underbrace{u_h-u^*|_h}_{e_{\rm sp}}.
$$

Measured algebraic error includes rounded input, reference construction,
and arithmetic effects; it is not pure exact-arithmetic iteration error.
Spatial error refers to the ideal discrete reference.
Reference samples/RHS are independently evaluated in long double and checked
against the public helpers. Historical tables remain in the
[research log](research-log.md).

## 5. Spatial Convergence

For a sampled error vector \(e\),

$$
\|e\|_{2,h}=\sqrt{h_xh_y}\left(\sum_{i,j}|e_{ij}|^2\right)^{1/2},
\quad \|e\|_\infty=\max_{i,j}|e_{ij}|,\quad
p_\ell=\frac{\log(E_{\ell-1}/E_\ell)}{\log(h_{\ell-1}/h_\ell)}.
$$

Both spacings halve, so the denominator is log(2). Squares use \(h=h_x=h_y\);
the rectangle uses \(h=\max(h_x,h_y)=h_y\). Orders are computed from total
errors separately in each norm. The first grid has no order (CSV value NA).

Taylor expansion independently gives

$$
\lambda_{mn}-\lambda_{mn,h}
=\frac{(m\pi/L_x)^4h_x^2+(n\pi/L_y)^4h_y^2}{12}
+O(h_x^4+h_y^4).
$$

A fixed positive gap and joint refinement predict second-order spatial error.
For the odd counts used here, the unit-amplitude mode satisfies
\(\|\phi_{11}\|_{2,h}=\sqrt{L_xL_y}/2,\ \|\phi_{11}\|_\infty=1\).
Single-mode spatial errors are \(d_{11,h}\sqrt{L_xL_y}/2\) and \(d_{11,h}\).
For the square mixture, sine orthogonality predicts
\(\|e_{\rm sp}\|_{2,h}=\frac12\sqrt{d_{11,h}^2+0.09d_{23,h}^2}\).
Its Linf reference is the maximum over actual nodes; the maximizing node can
change. Predicted orders use these references, never measured solver errors.

| Case | Interior grids | Spacings | Independent reference |
|---|---|---|---|
| Square single | 7x7,15x15,31x31,63x63 | h=1/8,1/16,1/32,1/64 m | Single-mode scale and norm identities |
| Square mixture | Same four square grids | Same h | Two modal scales, orthogonal L2, nodal Linf |
| Rectangle single | 7x5,15x11,31x23,63x47 | hx=1/8,1/16,1/32,1/64 m; hy=2hx | Geometry-specific modal scale and weight |

These are **recorded final-grid observations**, rounded from preserved CSVs,
not new measurements or cross-platform guarantees:

| Case | Iterations | True residual | Total L2 error | Total Linf error | Observed p L2 / Linf | Predicted p L2 / Linf |
|---|---|---|---|---|---|---|
| Square single | 1 | 1.335020e-13 | 1.057704e-4 | 2.115407e-4 | 2.000568 / 2.000568 | 2.000568 / 2.000568 |
| Square mixture | 3 | 9.159860e-12 | 2.501583e-4 | 6.082066e-4 | 2.003139 / 2.002857 | 2.003139 / 2.002857 |
| Rectangle single | 1 | 1.303567e-13 | 1.639120e-4 | 2.676671e-4 | 2.000692 / 2.000692 | 2.000692 / 2.000692 |

The rectangle's true_residual column uses the wide CSR product;
csr_true_residual records the fresh public matvec residual.
Single-mode RHS vectors are eigenvectors: one-step GMRES is expected and
provides no performance/scalability evidence. The mixture recorded 2,2,2,3
iterations; these convergence runs do not exercise restart.

## 6. Acceptance Criteria

These are existing policies, not new thresholds. Define
\(\varepsilon=\texttt{numeric\_limits<double>::epsilon()}\),
\(\eta=0.001\), and \(s_p=\log_2((1+\eta)/(1-\eta))\).

Stage 5A uses Lx=1 m and the following fixed inputs; nx and ny count interior
nodes. All cases start with amplitude 1+0.25i on the listed mode. The restarted
case adds modes (1,1) with amplitude 0.5-0.75i and (3,2) with amplitude -0.25+0.5i.

| Case suffix | nx,ny | Ly (m) | k (inverse metres) | Base mode |
|---|---|---|---|---|
| small | 3,2 | 1 | 3 | (1,2) |
| medium | 15,11 | 1 | 8 | (1,2) |
| altered | 7,5 | 0.75 | 3 | (2,1) |
| restarted | 5,4 | 1 | 3 | (1,2) |

| Experiment | Required conditions |
|---|---|
| Stage 5A | Zero start; tolerance=1e-12; restart=10 (single mode) or 2 (three-mode case); max_iterations=1000; CONVERGED; independent residual at most 1e-12; relative solution error at most \(2\kappa_2(10^{-12}+64\varepsilon)\); finite history of iterations+1 entries; global iteration limit; mixed-case iterations greater than restart. |
| All refinement studies | Zero start; tolerance=1e-11, restart=20, max_iterations=200; CONVERGED; \(0<\mathrm{iterations}\le200\); nonzero RHS; independent residual at most 1e-11; actual algebraic errors at most 0.001 times spatial errors in both norms; L2 budget below; residual at most 0.01 times relative spatial L2 error; strict decrease of total errors in both norms. |
| Square/rectangle single mode | Predicted p in [2,2.02]; observed p in [1.99,2.03]; each norm satisfies \(\lvert p-p_{\rm predicted}\rvert\le s_p\). |
| Square mixture | Each predicted and observed p in [1.99,2.06]; each norm satisfies \(\lvert p-p_{\rm predicted}\rvert\le s_p\); at least one grid takes more than one iteration. |

For the real symmetric positive-definite refinement operators, with ideal
minimum eigenvalue \(\delta_h\), tests use

$$
\|e_{\rm alg}\|_{2,h}\le B_2\le\eta\|e_{\rm sp}\|_{2,h},\qquad
B_2=\frac{\|r\|_{2,h}}{\delta_h}
       +64\varepsilon\kappa_2\|u_h\|_{2,h}.
$$

The inverse-gap term follows from the ideal SPD operator. The second term is
a conservative experimental roundoff allowance, not a rigorous interval bound
for every floating-point operation. The order perturbation allowance follows
from \((1-\eta)E_{\rm sp}\le E_{\rm total}\le(1+\eta)E_{\rm sp}\) on adjacent
grids. Broad p intervals additionally restrict the asymptotic regime.

Independent sample checks use absolute Linf limits 16 epsilon for single-mode
samples, 64 epsilon for mixture samples; RHS limits are 64 epsilon times the
sum of absolute continuous modal RHS amplitudes. Norm identities use
16 epsilon absolute error. These are checks at the recorded field scales,
not universal relative tolerances.

Thresholds are inclusive where coded as <= or interval endpoints; error decrease
is strictly <. Equality with the previous error fails.
NaN/Inf must not pass numerical acceptance. Drivers check finite complex
components before accumulation; residual/order predicates reject nonfinite
results, and the rectangle names 19 aggregate metrics in its finite-value loop.
Shared predicates assume finite, nonnegative norms and positive denominators;
they are not generic validators for arbitrary invalid synthetic inputs.
See the helper preconditions in tests/helmholtz/verification_checks.hpp.

Stage 5C.1 fixtures call the same test-only functions used by the square
drivers and applicable functions reused by the rectangle. They test valid
metrics, rejected statuses/residuals/errors/budgets/orders, nonfinite data,
inclusive boundaries and nextafter-outside values. A real two-mode solve
with max_iterations=1 must return NOT_CONVERGED after one iteration, with
independently verified residual above tolerance.
Typed AcceptanceFailure reasons are the primary contract. Tests also compare
the full condition field before "; measured[" and reject missing delimiters
or prefix-only matches. Unexpected exceptions and crashes do not count as
successful negative tests. These tests verify acceptance semantics, not
additional spatial convergence measurements.

## 7. Spectral Assumptions

For a uniform Dirichlet rectangle with at least one interior node per axis,
\(\sin z\ge2z/\pi\), \(0\le z\le\pi/2\), gives

$$
\lambda_{11,h}\ge \frac4{L_x^2}+\frac4{L_y^2},\qquad
\delta_h=\lambda_{11,h}-k^2.
$$

For k=1 inverse metre in the refinement studies the full spectrum is positive:
\(\delta_h\ge7\,\mathrm{m}^{-2}\) on the square and
\(\delta_h\ge(43/9)\,\mathrm{m}^{-2}\) on the rectangle.
Every discrete mode lies above the first eigenvalue, so these bounds cover
the full spectrum, not just the manufactured mode.
The largest eigenvalue occurs at (nx,ny); hence
\(\kappa_2=(\lambda_{n_xn_y,h}-k^2)/\delta_h\).
The rectangle uses its own bound, not the square-specific predicate.

Stage 5A instead enumerates all discrete modal eigenvalues, requiring
\(\min|\lambda_{mn,h}-k^2|>1\,\mathrm{m}^{-2}\) and \(\kappa_2<1000\).
Its k values are 3 or 8 inverse metres, not uniformly k=1; positive
definiteness is not assumed for every integration case.
These checks do not replace continuous uniqueness checks for a future PDE
experiment. Generic ill-conditioned GMRES tests are not a near-resonant
Helmholtz study.

## 8. Reproducibility

Use the baseline above for code and historical CSV evidence; this page is a
documentation-only addition. Stage-specific provenance remains in the
[research log](research-log.md). Relevant code checkpoints are:

| Stage | Commit |
|---|---|
| 5A integration | 69c27752fc022a04d200fa4f176dc9011db65b59 |
| 5B square convergence | 9c22d972af32ff13c43cb9d5d2b62f63d5bcb3ce |
| 5C.1 acceptance semantics | b99f45ae8e056efe2582ddf6186db48b3bf4aafe |
| 5C.2 rectangular convergence | 58fcf252f97e69b60b6bf1efe2fd0ece9286c8a6 |

The recorded environment is Ubuntu/WSL2
x86-64, GCC 13.3.0, CMake 3.28.3, C++20 Debug with
-g -Wall -Wextra -Wpedantic. Production values are complex double;
references/norms use long double where specified. Its precision is platform
dependent. No random matrices or external C++ libraries are required.

From the repository root with existing prerequisites installed:

~~~bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure
ctest --test-dir build -R '^wavecore.integration\.' -V
ctest --test-dir build -R '^wavecore.convergence' -V
ctest --test-dir build -R '^wavecore.verification.acceptance$' -V
.venv/bin/python -m mkdocs build --strict
git diff --check
~~~

The three preserved CSV files under experiments/helmholtz_2d/results/ are:

- stage5b_grid_convergence.csv — square single mode;
- stage5b_multimode_convergence.csv — square mixture;
- stage5c_rectangular_convergence.csv — rectangular single mode.

They are version-controlled historical evidence. Executables print a metadata
line followed by CSV and do not overwrite these files. Compare repeat output
without modifying data:

~~~bash
.venv/bin/python - <<'PY'
from pathlib import Path
import subprocess

cases = (
    ("wavecore_convergence_test", "stage5b_grid_convergence.csv"),
    ("wavecore_multimode_convergence_test", "stage5b_multimode_convergence.csv"),
    ("wavecore_rectangular_convergence_test", "stage5c_rectangular_convergence.csv"),
)
for executable, filename in cases:
    runs = [subprocess.check_output([str(Path("build") / executable)]) for _ in range(2)]
    expected = (Path("experiments/helmholtz_2d/results") / filename).read_bytes()
    assert runs[0] == runs[1], executable
    assert runs[0].split(b"\n", 1)[1] == expected, filename
    print(filename, "byte-identical")
PY
~~~

Byte identity is a same-environment regression check. Other compilers, standard
libraries, formats, flags or CPUs can change final digits or iteration counts;
cross-platform agreement has not been verified. Investigate mismatches rather
than silently replacing CSVs or weakening thresholds. Record new environment
metadata and compare mathematical criteria independently.
See [developer commands](developer-commands.md) for prerequisites and separate
local CI commands. GitHub Actions are not used.

## 9. Verification Traceability Matrix

Evidence paths are repository-relative. CSV names refer to the three files
above. Complete per-test details remain in the [coverage matrix](test-coverage.md).

| Mathematical claim / property | Independent reference | Test / CTest name | CSV or evidence | Acceptance criterion | Limitation |
|---|---|---|---|---|---|
| Interior geometry and x-fastest mapping | Explicit coordinates and inverse mapping | wavecore.grid | tests/helmholtz/test_grid.cpp | Coordinate tolerance, bijection, typed invalid-input/overflow failures | Finite set of grids |
| Complex sparse product and CSR invariants | Hand-computed complex products | wavecore.csr_matrix | tests/helmholtz/test_csr_matrix.cpp | Expected products and structural rejection assertions | No broad arithmetic-scale survey |
| Five-point stencil, BC omission and no row wrapping | Explicit small-row stencils and dense transpose | wavecore.assembly | tests/helmholtz/test_assembly.cpp | Small-matrix comparisons; nnz=1,12,20,105; symmetry and sorted columns | Constant coefficients, uniform grids |
| Continuous/discrete modal RHS | Sine formulas and CSR modal identity | wavecore.assembly; all three convergence tests | Assembly assertions; three CSVs | Modal matvec identity and continuous RHS limits in Section 6 | Selected modes only |
| Complex Arnoldi/GMRES diagnostics | Known solutions, wide residual, one-pass control, basis Gram matrix | wavecore.gmres | tests/helmholtz/test_gmres.cpp; campaign in research log | Case-specific error/residual/orthogonality/status assertions | Not a universal stability proof |
| Algebraic integration and restart | Full modal spectrum, sampled solution, fresh CSR residual | wavecore.integration.small / medium / altered / restarted | tests/helmholtz/test_integration.cpp; Stage 5A log | Stage 5A criteria in Section 6; restarted iterations>2 | Four small systems; no spatial claim |
| Square single-mode spatial order | Discrete scaling, sine norm identities | wavecore.convergence | stage5b_grid_convergence.csv | Single-mode p intervals, prediction agreement, algebraic separation | Eigenvector RHS |
| Square mixture spatial order | Two scales, orthogonal L2 and nodal Linf | wavecore.convergence.multimode | stage5b_multimode_convergence.csv | Mixture p intervals, norm-specific predictions; multiple iterations | Two modes; restart not exercised |
| Rectangular geometry, norms and order | Physical-coordinate formulas and wide CSR product | wavecore.convergence.rectangular | stage5c_rectangular_convergence.csv | Joint halving, geometry/gap checks, both residuals, single-mode p criteria | One aspect ratio; fixed hx/hy |
| Solver/spatial error separation | Ideal discrete u_h and residual/gap scale | All three convergence tests | Algebraic/spatial/budget CSV columns | Both algebraic norms and L2 budget below 0.1% spatial errors | Roundoff budget not universal |
| Acceptance semantics and finite/boundary rejection | Synthetic fixtures, nextafter, real limited solve | wavecore.verification.acceptance | tests/helmholtz/test_verification_checks.cpp | Correct typed reason and full condition; no arbitrary exception accepted | Documented preconditions; combinations not exhaustive |
| Build identity (non-numerical) | CMake expected metadata | wavecore.library_metadata | tests/helmholtz/test_library.cpp | Exact project/version checks | No PDE claim |

## 10. Known Limitations and Future Work

Verified scope is a few small algebraic systems, two square manufactured
fields, and one rectangular mode at fixed spacing ratio. References share the
ideal uniform-grid model; modeling errors are outside these tests.
Arbitrary forcing, irregular geometry, unequal-rate refinement, heterogeneous
media, PML and application-level solve I/O are not covered.

No cross-platform verification, performance/scaling evidence, or systematic
Helmholtz resonance/conditioning sweep is established. Small standalone GMRES
diagnostics do not fill those PDE-level gaps. No preconditioner is implemented.

Future work includes Stage 6 spectral-distance sweeps with residual and forward
error measured separately, then justified preconditioning and matched
performance experiments. PETSc benchmarking is a possible future independent
comparison, not an implemented dependency or completed validation. It requires
approved dependency integration, matched matrices/RHS and stopping criteria,
independent residual checks, and recorded versions, hardware, flags, setup
cost, runtime and memory. See the [roadmap](research-roadmap.md).

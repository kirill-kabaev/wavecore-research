# C++ API

The public numerical API is declared in `include/wavecore/helmholtz/` and uses
the `wavecore::helmholtz` namespace. The library is C++20 and has no external
C++ dependencies.

## `Grid2D`

`Grid2D` represents a uniform rectangular grid containing only interior
unknowns. Its constructor is:

```cpp
Grid2D(std::size_t nx, std::size_t ny, double lx, double ly);
```

`nx` and `ny` are the numbers of interior nodes, while `lx` and `ly` are the
positive domain lengths. The spacings are

```text
hx = lx / (nx + 1),    hy = ly / (ny + 1).
```

The public query and mapping functions are:

```cpp
std::size_t nx() const noexcept;
std::size_t ny() const noexcept;
std::size_t unknown_count() const noexcept;
double lx() const noexcept;
double ly() const noexcept;
double hx() const noexcept;
double hy() const noexcept;
double x(std::size_t i) const;
double y(std::size_t j) const;
std::size_t linear_index(std::size_t i, std::size_t j) const;
Grid2D::index_pair grid_indices(std::size_t p) const;
```

The zero-based API indices correspond to the mathematical interior nodes
`x_(i+1)` and `y_(j+1)`. Coordinates are `(i + 1) * hx` and `(j + 1) * hy`.

### X-fastest indexing

Linear storage uses

```text
p = j * nx + i.
```

The x index `i` therefore changes in consecutive memory locations. The inverse
mapping is `i = p % nx`, `j = p / nx`. This ordering is also used by operator
assembly and manufactured vectors.

Zero grid dimensions, nonfinite or nonpositive lengths, and invalid computed
spacings cause `std::invalid_argument`. Size calculations that cannot be
represented by `std::size_t` cause `std::overflow_error`. Coordinate and index
queries outside the grid cause `std::out_of_range`.

## `CsrMatrix`

`CsrMatrix` is an immutable compressed sparse row matrix. Its public aliases
are:

```cpp
using size_type = std::size_t;
using value_type = std::complex<double>;
using index_array = std::vector<size_type>;
using value_array = std::vector<value_type>;
```

The constructor receives the complete CSR representation:

```cpp
CsrMatrix(size_type rows,
          size_type cols,
          index_array row_ptr,
          index_array col_idx,
          value_array values);
```

Complex entries are stored as `std::complex<double>`; both their real and
imaginary parts must be finite. The three arrays have the usual CSR meaning:

- `row_ptr` has `rows + 1` entries and identifies the half-open entry range
  `[row_ptr[r], row_ptr[r + 1])` for row `r`;
- `col_idx[q]` is the column of stored entry `q`;
- `values[q]` is its complex value.

The representation requires `row_ptr[0] == 0`, `row_ptr.back() == nnz`,
monotonically nondecreasing row pointers, in-range column indices, and strictly
increasing column indices within every row. Duplicate column entries are
therefore rejected. Empty matrices and empty rows are supported when these
invariants hold.

The read-only API is:

```cpp
size_type rows() const noexcept;
size_type cols() const noexcept;
size_type nnz() const noexcept;
const index_array& row_ptr() const noexcept;
const index_array& col_idx() const noexcept;
const value_array& values() const noexcept;
value_array multiply(const value_array& vector) const;
```

`multiply` returns a newly allocated vector of length `rows()`. Its input must
have length `cols()` or it throws `std::invalid_argument`. Invalid CSR input
throws `std::invalid_argument`; an unrepresentable `rows + 1` size throws
`std::overflow_error`.

## 2D Helmholtz FDFD assembly

```cpp
CsrMatrix assemble_helmholtz_operator(const Grid2D& grid,
                                      double wavenumber);
```

This function builds the real-valued five-point stencil inside the complex CSR
container for

```text
A = -Delta_h - k^2 I
```

with homogeneous Dirichlet boundary conditions. Its coefficients are

```text
diagonal:       2 / hx^2 + 2 / hy^2 - k^2
x-neighbours:  -1 / hx^2
y-neighbours:  -1 / hy^2.
```

Boundary values do not appear as unknowns, so missing neighbours are omitted.
The resulting square matrix uses the same x-fastest ordering as `Grid2D`, and
its column indices are strictly increasing within every row.

`wavenumber` must be finite and nonnegative. Invalid values cause
`std::invalid_argument`; nonrepresentable squared wavenumbers, stencil
coefficients, or CSR sizes cause `std::overflow_error`. An internal mismatch
between the predicted and assembled nonzero counts causes `std::logic_error`.

## Manufactured modes and right-hand sides

The modal and manufactured-data functions are:

```cpp
double continuous_mode_eigenvalue(const Grid2D& grid,
                                  std::size_t mode_x,
                                  std::size_t mode_y);

double discrete_mode_eigenvalue(const Grid2D& grid,
                                std::size_t mode_x,
                                std::size_t mode_y);

CsrMatrix::value_array manufactured_solution(
    const Grid2D& grid,
    std::size_t mode_x,
    std::size_t mode_y,
    CsrMatrix::value_type amplitude = {1.0, 0.0});

CsrMatrix::value_array continuous_manufactured_rhs(
    const Grid2D& grid,
    std::size_t mode_x,
    std::size_t mode_y,
    double wavenumber,
    CsrMatrix::value_type amplitude = {1.0, 0.0});

CsrMatrix::value_array algebraic_manufactured_rhs(
    const Grid2D& grid,
    std::size_t mode_x,
    std::size_t mode_y,
    double wavenumber,
    CsrMatrix::value_type amplitude = {1.0, 0.0});
```

For positive integer modes `m = mode_x`, `n = mode_y`, the sampled solution is

```text
u(i,j) = amplitude
       * sin((i + 1) m pi / (nx + 1))
       * sin((j + 1) n pi / (ny + 1)).
```

The continuous eigenvalue is
`lambda_mn = (m*pi/lx)^2 + (n*pi/ly)^2`. The discrete eigenvalue uses the
finite-difference spectrum,
`lambda_mn,h = 4 sin^2(m*pi/(2(nx+1))) / hx^2`
`+ 4 sin^2(n*pi/(2(ny+1))) / hy^2`.

The two right-hand sides serve different checks:

- `continuous_manufactured_rhs` returns `(lambda_mn - k^2) * u`. It represents
  the continuous PDE sampled on the grid and therefore includes discretization
  error when paired with the discrete operator.
- `algebraic_manufactured_rhs` returns `(lambda_mn,h - k^2) * u`. It is
  algebraically consistent with the assembled finite-difference operator, up
  to floating-point roundoff.

Modes must be positive and no larger than the corresponding number of interior
nodes. The amplitude must have finite real and imaginary parts, and the
wavenumber must be finite and nonnegative. Invalid inputs cause
`std::invalid_argument`; nonfinite computed values cause `std::overflow_error`.

## Restarted complex GMRES

The solver entry point and option types are:

```cpp
struct GmresOptions {
    double tolerance = 1.0e-10;
    std::size_t restart = 30;
    std::size_t max_iterations = 1000;
};

enum class GmresStatus { CONVERGED, NOT_CONVERGED, BREAKDOWN };

GmresResult gmres(
    const CsrMatrix& matrix,
    const CsrMatrix::value_array& rhs,
    const GmresOptions& options = {},
    const CsrMatrix::value_array& initial_guess = {});
```

The implementation is sequential and unpreconditioned. It uses the Hermitian
inner product `sum(conj(v[i]) * w[i])`, two modified Gram--Schmidt passes, and
complex Givens rotations. An empty `initial_guess` selects the zero vector.

- `tolerance` is the positive finite stopping threshold applied to the true
  residual measure.
- `restart` is the maximum Krylov dimension per cycle. It must be positive and
  is capped at the system dimension.
- `max_iterations` is a positive global limit across all restart cycles.

The matrix must be square; `rhs` and a nonempty initial guess must match its
dimension. Options and input vector values are validated. Invalid dimensions,
options, or nonfinite inputs cause `std::invalid_argument`. Nonfinite arithmetic
causes `std::overflow_error`.

### Result and statuses

```cpp
struct GmresResult {
    CsrMatrix::value_array solution;
    GmresStatus status;
    std::size_t iterations;
    double final_residual;
    std::vector<double> residual_history;
};
```

- `CONVERGED` means a recomputed true residual satisfies `tolerance`.
- `NOT_CONVERGED` means the global iteration limit was reached without meeting
  the tolerance. Restart-induced stagnation has no separate status.
- `BREAKDOWN` means that no further Krylov direction was available, or the
  triangular least-squares solve was singular, before true convergence was
  confirmed.

For nonzero `rhs`, `final_residual` is `||b - A*x||_2 / ||b||_2`. If `rhs` is
exactly zero, it is the absolute norm `||b - A*x||_2`. `residual_history`
contains the true residual at the initial guess followed by one entry per
Arnoldi iteration. The least-squares residual estimate does not by itself
produce `CONVERGED`; the implementation explicitly recomputes `b - A*x`.

A small residual does not guarantee a small forward solution error for an
ill-conditioned matrix. Near-breakdown, finite-precision cancellation,
underflow, overflow, and restart stagnation remain numerical limitations. The
API makes no universal stability or scale-invariance guarantee.

## Minimal composition example

The following code uses only public declarations. It manually composes grid
construction, operator assembly, manufactured algebraic data, and GMRES; it
does not represent a dedicated end-to-end application workflow.

```cpp
#include "wavecore/helmholtz/assembly.hpp"
#include "wavecore/helmholtz/gmres.hpp"

using namespace wavecore::helmholtz;

Grid2D grid{31, 23, 1.0, 0.75};
const double wavenumber = 8.0;

const CsrMatrix matrix = assemble_helmholtz_operator(grid, wavenumber);
const auto rhs = algebraic_manufactured_rhs(grid, 1, 2, wavenumber);

GmresOptions options;
options.tolerance = 1.0e-10;
options.restart = 30;
options.max_iterations = 500;

const GmresResult result = gmres(matrix, rhs, options);
if (result.status == GmresStatus::CONVERGED) {
    // result.solution contains the computed complex-valued vector.
}
```

The caller remains responsible for interpreting convergence, conditioning,
discretization error, physical units, and resonance proximity.

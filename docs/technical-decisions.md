# Technical Decisions

This page explains choices evidenced by the current Stage 5A/5B sources and
[experiment records](research-log.md). It documents the current rationale;
it does not reconstruct undocumented historical discussions or approval dates.

## Continuous manufactured forcing for spatial verification

The model is the homogeneous, isotropic, lossless TM_z Helmholtz problem
with time convention exp(-i omega t), zero Dirichlet/PEC walls, and
A_h=-Delta_h-k squared I. The experiments use a unit square in metres and
k=1 inverse metre. Grid dimensions count interior nodes, h=1/(n+1), with
x-fastest storage p=j*nx+i.

For a smooth zero-boundary sine mode phi_mn, continuous forcing is
f=(lambda_mn-k squared) phi_mn, where
lambda_mn=pi squared (m squared/Lx squared + n squared/Ly squared).
Stage 5B uses the existing `continuous_manufactured_rhs` helper and checks
its output against independent coordinate-based long-double formulas.

Algebraic forcing b=A_h*u_exact makes the nodal exact field solve the discrete
system by construction. That is appropriate for Stage 5A integration, but
removes the spatial error that Stage 5B needs to measure. Keeping continuous
forcing exposes the difference between continuous and discrete eigenvalues.

## Separate numerical/algebraic and discretization errors

The reference quantities are u*, the continuous field sampled at interior nodes,
and u_h, the exact-arithmetic solution of the discrete stencil with continuous
forcing. The error decomposition is

$$
u_{\mathrm{computed}}-u^*|_h
=(u_{\mathrm{computed}}-u_h)+(u_h-u^*|_h).
$$

The first term includes solver error and floating-point effects in sampled
inputs and arithmetic; it is not a pure exact-arithmetic iteration error. The
second is the spatial discretization error. Both are measured in absolute
weighted L2(h)=sqrt(hx hy sum |e| squared) and Linf=max |e|.

Both tests require actual numerical/algebraic errors <=0.1% of spatial errors
in both norms. They also check a conservative residual/rounding-aware L2
budget. This prevents a small residual alone from being interpreted as adequate
forward accuracy. The residual is recomputed using a fresh public CSR matvec
and wide norm accumulation, independently of GMRES residual bookkeeping;
the matvec kernel itself is shared with the production solver.

## Analytic discrete references

On the uniform Dirichlet grid, sampled sine modes are eigenvectors of the
five-point stencil:

$$
\lambda_{mn,h}
=\frac4{h_x^2}\sin^2\!\left(\frac{m\pi h_x}{2L_x}\right)
+\frac4{h_y^2}\sin^2\!\left(\frac{n\pi h_y}{2L_y}\right).
$$

For amplitudes a_mn and continuous modal forcing,

$$
u_h=\sum_{m,n}a_{mn}
\frac{\lambda_{mn}-k^2}{\lambda_{mn,h}-k^2}\phi_{mn}|_h.
$$

The tests compute these references through independent long-double trigonometric
formulas, not by invoking the production discrete-eigenvalue helper or another
linear solver. This avoids introducing a second numerical solver and allows
direct separation of errors. Long double remains finite-precision arithmetic;
these checks are not a proof for all grids or platforms.

## Complementary single-mode and multimode experiments

The single-mode field phi_11 gives a simple closed-form error factor and known
discrete norms. Its RHS is a discrete eigenvector, so one GMRES iteration is
expected in exact arithmetic. It is a clear spatial test but a weak integration
challenge for the Krylov solver.

The separate field phi_11 - 0.3 phi_23 excites two distinct eigenvalues. It
retains an exact discrete reference while requiring more than one iteration.
The test asserts multiple iterations on at least one grid, not exactly two:
rounding can excite additional components. Its Linf maximum can move between
nodes, so finite-grid orders need not approach two monotonically.

Both experiments halve h exactly using n=7,15,31,63. Acceptance intervals are
tied to independently computed finite-grid modal predictions and the 0.1%
algebraic-error budget. The resulting allowable order perturbation is
log2(1.001/0.999). This is more specific than assuming every measured order must
equal two exactly. The single-mode and multimode CSV files remain separate.

## Scope of the spectral checks

For k=1 on the unit square the operators are positive definite. The inequality
sin(z)>=2z/pi gives lambda_11,h>=8, hence the spectral gap is at least 7 inverse
square metres. Extreme sine-mode eigenvalues give kappa_2 for these real
symmetric matrices. The tests use this gap to relate a residual to an L2
forward-error scale and to check adequate separation from spatial error.

These spectral formulas depend on the rectangular, uniform, constant-coefficient
Dirichlet model. They do not validate irregular grids, heterogeneous media,
radiation boundaries, PML, or general nonnormal complex operators. The
64 epsilon kappa_2 rounding allowance is an experiment-specific conservative
budget, not a universal numerical-stability theorem. Near-resonance studies,
preconditioning, and scalability remain separate roadmap work.

## Verification and diagnostics

The [coverage matrix](test-coverage.md) lists the registered tests. The
[research log](research-log.md#stage-5b-grid-convergence-study) and
[multimode record](research-log.md#stage-5b-hardening-multimode-convergence)
preserve parameters, environment, full CSV paths, results, and reproduction
commands. Failure diagnostics report n, h, the condition, measured quantities,
and allowed limits without altering the scientific acceptance criteria.

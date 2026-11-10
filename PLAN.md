# PLAN: Add constant-magnitude constraint and Tikhonov smoothing to the existing Split-Bregman TV reconstruction

## 0. Goal

The repository already reconstructs the magnetization vector field **m** (3 components on an N-voxel grid) from multi-axis laminography XMCD data. It uses Split-Bregman (ADMM) with a TV regularizer. This plan extends that solver with:

1. **A Tikhonov gradient penalty** `(β/2)‖∇m‖²`.
2. **A constant-magnitude constraint** |m(r)| = Ms on the magnetic support Ω, enforced through a second split `z = m`. The constraint set is a sphere at each voxel, which is nonconvex, so the solver converges to a stationary point only.
3. Penalty continuation, optional estimation of Ms, and diagnostics.

**Hard requirement:** with `beta = 0` and the magnitude constraint off (`use_magnitude = False`, equivalently ρ = 0), the solver must reproduce the current TV results, within floating-point tolerance of the existing output. Every new feature is opt-in.

---

## 1. Reference math (implement exactly this)

### Notation
- `m ∈ R^{3×N}`: magnetization.
- `∇`: the gradient applied separately to each component. `(∇m)(r) ∈ R^{3×D}`, where D is the number of spatial dimensions.
- `A`: the forward operator. `AᵀA` is the existing normal operator (possibly Toeplitz/FFT based).
- `Ω`: the support mask. `g`: the data.
- Bregman sign convention (Goldstein–Osher): `b ← b + (∇m − d)`.

### Objective
```
min_m  ½‖A m − g‖² + λ Σ_r ‖(∇m)(r)‖_F + (β/2)‖∇m‖²   s.t.  m ∈ S
S: |m(r)| = Ms on Ω,  m(r) = 0 off Ω
```

### Split
`d = ∇m` (penalty μ, Bregman variable `b_d`) and `z = m` (penalty ρ, Bregman variable `b_z`).

### Iteration k → k+1
**(b) m-step: linear solve**
```
(AᵀA + (μ+β) ∇ᵀ∇ + ρ I) m^{k+1} = Aᵀg + μ ∇ᵀ(d^k − b_d^k) + ρ (z^k − b_z^k)
```
- β appears **only in the system matrix**, never on the right-hand side.
- ρ adds a diagonal term. The matrix is SPD for ρ > 0, even when A has a null space.

**(c) d-step: TV shrinkage.** Keep whichever variant the code already uses. The vector-TV form is:
```
w(r) = (∇m^{k+1})(r) + b_d^k(r)
d^{k+1}(r) = max(‖w(r)‖_F − λ/μ, 0) · w(r)/‖w(r)‖_F      (0 if w(r) = 0)
```

**(d) z-step: projection onto the sphere, independent per voxel.** Let `v = m^{k+1} + b_z^k`.
```
z(r) = Ms · v(r)/|v(r)|     on Ω, |v(r)| > eps
z(r) = z^k(r)               on Ω, |v(r)| ≤ eps   (projection not unique: keep previous)
z(r) = 0                    off Ω
```
The d-step and z-step both depend only on m^{k+1}. They are independent of each other and can run in either order.

**(e) Bregman updates**
```
b_d^{k+1} = b_d^k + ∇m^{k+1} − d^{k+1}
b_z^{k+1} = b_z^k + m^{k+1} − z^{k+1}
```

### Penalty continuation
When ρ → ρ', rescale `b_z ← (ρ/ρ') b_z`. Apply the same rule to μ and `b_d` if μ is ever changed. This keeps the unscaled multiplier fixed.

### Optional Ms estimation
After computing v and before projecting:
```
Ms = mean_{r∈Ω} |v(r)|
```
This is the closed-form least-squares estimate of a scalar Ms.

### Residuals (normalize by ‖m‖ or ‖z‖)
```
primal TV:   ‖∇m − d‖          dual TV:   μ ‖∇ᵀ(d^{k+1} − d^k)‖
primal mag:  ‖m − z‖           dual mag:  ρ ‖z^{k+1} − z^k‖
```

### Output
Return **z** as the reconstruction, since it satisfies |z| = Ms exactly on Ω. Also keep **m** and the final `b_z` for diagnostics.

### Convergence facts (do not "fix" behavior that matches these)
- The sphere constraint set is nonconvex, so the standard ADMM convergence proof does not apply.
- The available guarantee is convergence to a **stationary point**, and only for ρ large enough (Wang, Yin & Zeng 2019). There is no global-optimality guarantee.
- Expected behavior includes dependence on the initialization and non-monotone residuals at fixed ρ. Increasing ρ during the solve (continuation) and warm-starting from the unconstrained TV solution are the mitigations.
- Without ρ, the TV split alone is convex and converges globally (unchanged from the current code).

---

## 2. Ground rules for the implementer

- **Phase 0 (discovery) comes before any edits.** Write the findings into a `## Findings` section appended to this file.
- Do not change the existing discretization, boundary conditions, forward model, data pipeline, or default parameters. If something looks wrong (for example central differences, or an adjoint that isn't exact), **record it in Findings and ask**. Do not silently fix it.
- Follow existing conventions: array layout, backend (numpy/cupy/torch), dtype, config/parameter system, logging, naming.
- Add no new dependencies.
- Make small commits, one per phase, each with passing tests.

---

## Phase 0: Discovery (read-only)

Locate the following and record each answer in Findings, with file paths and function names:

1. **Solver entry point**: the TV Split-Bregman/ADMM loop, its parameters (λ, μ, iteration count, tolerances) and what it returns.
2. **Forward/adjoint and normal operator**: how `AᵀA` is applied (direct NUFFT, Toeplitz/FFT, explicit matrix), and where `Aᵀg` is computed.
3. **Gradient and divergence**: the type of finite differences, the boundary condition, and whether `div = −∇ᵀ` exactly. Confirm by running a dot-product test `⟨∇x, P⟩ = ⟨x, ∇ᵀP⟩` on random arrays.
4. **m-step solver**: CG (tolerance and max iterations), FFT/DCT direct solve, or Gauss–Seidel; whether it uses a preconditioner; how the system operator is assembled.
5. **Shrink variant**: vector (Frobenius over components and directions), per-component isotropic, or anisotropic.
6. **Bregman variable convention**: the sign and scaling of `b_d`, and whether μ is constant or continued.
7. **Array layout and dtype**: for example `(3, Nz, Ny, Nx)`, real or complex, CPU or GPU, single or double precision.
8. **Support mask**: whether a mask Ω exists and how it is applied. If none exists, note that one must be added as an input, defaulting to all True.
9. **Config/CLI**: how solver parameters are passed and stored.
10. **Tests**: the existing test framework, fixtures, and synthetic phantoms or forward simulations that can be reused.
11. **Memory**: z and `b_z` each add one array the size of m. Estimate the extra memory at production size and note whether it fits.

**Exit criterion:** Findings is complete. The adjoint test passes, or its failure is reported to the user.

---

## Phase 1: Tikhonov term (β)

1. Add parameter `beta` (default `0.0`).
2. Add `β ∇ᵀ∇` to the m-step system operator:
   - **CG:** add `beta * gradT(grad(x))`, or fold it in as `(mu + beta)` multiplying `gradT(grad(x))`.
   - **FFT/DCT direct solve:** add `β|D̂(k)|²` to the denominator, where `|D̂(k)|²` is the symbol already used for μ.
   - **Preconditioner:** update it the same way if it includes the μ term.
3. Leave the right-hand side unchanged.
4. Tests:
   - `beta=0` reproduces the current output (regression).
   - The objective with β > 0 is evaluated correctly.
   - The m-step operator stays symmetric (dot-product test).

---

## Phase 2: Magnitude split (z, b_z)

1. Add parameters:
   - `use_magnitude` (bool, default `False`)
   - `rho` (float)
   - `Ms` (scalar, or per-voxel array if trivial to support)
   - `support` (bool mask over the grid; default all True)
2. **State:** add `z` and `b_z`, each shaped like `m`. If `use_magnitude` is False, allocate neither and skip every new term. This guarantees the regression requirement.
3. **Initialization:**
   - Accept an optional `m0` so the solver can be warm-started (intended use: the unconstrained TV result).
   - `z⁰ = P_S(m⁰)` and `b_z⁰ = 0`.
   - Where `|m⁰(r)| ≤ eps` on Ω, use a fixed default direction (for example +x) or an optional caller-supplied `z0`.
4. **m-step:** add `rho * x` to the system operator and `rho * (z − b_z)` to the right-hand side. With an FFT direct solve, add `ρ` to the denominator. ρI is diagonal in Fourier space, so the existing diagonalization still holds.
5. **z-step:** implement one pure function, `project_sphere(v, Ms, support, z_prev, eps)`, vectorized over voxels. `eps` guards against division by zero.
6. **Bregman update** for `b_z` as in §1(e).
7. **Return values:** z (the reconstruction), m, `b_z`, and the history.

---

## Phase 3: Continuation and Ms estimation

1. **ρ schedule:** optional `rho_growth` (default 1.0) and `rho_max`. Each time ρ changes, rescale `b_z ← (ρ_old/ρ_new) b_z`. If the m-step uses a direct FFT solve or a precomputed preconditioner, refresh it when ρ changes.
2. **Only monotone increases of ρ.** Do not implement residual-balancing ρ adaptation that can decrease ρ; its convergence theory covers convex problems only.
3. **Optional `estimate_Ms`:** `Ms ← mean_{Ω}|m + b_z|` before projecting. Log its value every iteration.

---

## Phase 4: Diagnostics and logging

Log at each logging interval:
- the four residuals in §1
- the individual objective terms: data misfit, TV, Tikhonov
- ρ, and Ms if estimated

Save at the end:
- the `|m|/Ms` map and its histogram over Ω, where m is the unprojected ADMM variable. Voxels far from 1 are where the constraint fights the data: domain walls, vortex cores, or components the operator cannot see.
- the multiplier map `ρ|b_z|` over Ω, which shows where the data push |m| away from Ms
- the per-iteration history (a dict or list) returned by the solver

---

## Phase 5: Tests

### Unit tests
- **Projection:**
  - Idempotent: `P(P(v)) = P(v)`.
  - Output is zero off the support.
  - `|z| = Ms` on Ω, to within 1e-12 relative error in double precision.
  - Direction is preserved: `z(r) ∥ v(r)` with a positive scale factor.
  - The zero-vector case keeps the previous z.
- **Adjoint test** for `grad`/`gradT`, to about 1e-12 relative error in double precision.
- **m-step operator** is symmetric and positive: `⟨x, Hy⟩ = ⟨Hx, y⟩` and `⟨x, Hx⟩ > 0`.
- **ρ rescaling** leaves `ρ·b_z` unchanged.

### Regression test
- With `use_magnitude=False` and `beta=0`, output matches the current code on the existing test case.

### Synthetic end-to-end tests
Use a small grid, for example 32³ or a 2-D slice.
1. **Denoising (A = I):** a unit-magnitude domain phantom plus noise.
   - Error falls below the noisy-input error.
   - `max_Ω ||z| − Ms| < 1e-6`.
   - Primal and dual magnitude residuals reach 1e-4 with ρ continuation.
2. **Initialization sensitivity (report, do not assert equality):** run from three random starts plus the unconstrained-TV warm start. Report the final objective, the error against ground truth, and the spread across starts. This documents the nonconvexity. The test passes if the warm start reaches an objective value no worse than the best random start.
3. **Production forward operator on a small phantom:** compare TV only and TV + magnitude constraint against ground truth. Report error and iteration counts.

---

## Phase 6: Recommended usage (put in a docstring or README section)

1. Run TV (plus β if used) with no constraint, to get a baseline and `m0`.
2. Run with `use_magnitude=True`, warm-started from `m0`, with ρ continuation.
3. Inspect the `|m|/Ms` and `ρ|b_z|` maps to see where the constraint is doing the most work.

Parameter starting points. **These are heuristics, not established results.** Document them as such.
- `rho`: start at the same order as μ. If residuals stall or oscillate, use `rho_growth ≈ 1.01–1.05` up to `rho_max ≈ 10–50 × rho0`.
- `beta`: small relative to μ. Pick it with the discrepancy principle: the data misfit should be about `½ · (number of measurements) · σ²`, where σ is the noise standard deviation. Don't tune it by eye.
- `Ms`: in XMCD units it is usually unknown. Use `estimate_Ms`, or a high percentile of `|m0|` on Ω, and report which was used.

---

## Out of scope
- A relaxed |m| ≤ Ms (ball) constraint.
- Changing the TV variant, the discretization, the boundary conditions, or the forward model.
- Tuning λ or μ, or retuning production parameters.
- GPU porting or other performance work beyond what the new terms require.
- Changes to the data pipeline or the ptychography preprocessing.

## Deliverables checklist
- [ ] `## Findings` section appended to this file (Phase 0)
- [ ] β term with passing regression test (Phase 1)
- [ ] Sphere projection plus the z/b_z split, with unit tests (Phase 2)
- [ ] ρ continuation with dual rescaling and optional Ms estimation (Phase 3)
- [ ] Residual and objective logging plus magnitude and multiplier maps (Phase 4)
- [ ] Synthetic end-to-end tests passing (Phase 5)
- [ ] Usage documentation with heuristics clearly labeled (Phase 6)
- [ ] Short summary for the user: files changed, new parameters, test results, open questions from Findings

---

## Findings

Phase 0 results (2026-10-09). Paths are relative to the repo root.

1. **Solver entry point.** CPU `opt::split_bregman`, `src/bregman.cpp:38-137`; GPU `gpu::opt::split_bregman`, `src/gpu/bregman.cu`. Parameters: lambda, mu (constant), outer/inner iterations, tol, xtol, support. Returns x. The constrained/Tikhonov solver is a separate function, `opt::regularized_admm` (`src/regularized_admm.cpp`, `src/gpu/regularized_admm.cu`); `split_bregman` is untouched.
2. **Normal operator.** `A^T A` is a `Function<T>` built in `src/mbir.cpp` / `src/gpu/mbir.cu` from Toeplitz/FFT kernels (`cpu::sysmat`). `A^T g` (`yT`) is the back-projected data computed before the solver call.
3. **Gradient / divergence.** `src/finitediff.cpp`, `src/gpu/finitediff.cu`. HEAD used central differences with an inexact adjoint. The working tree uses forward differences with Neumann (zero difference at the last index); `divergence` is the exact adjoint (`test_regularized_admm`: dot-product test, 1e-9; GPU test likewise). **Deviation from the "do not change the discretization" rule; approved by the user (2026-10-09).** Regression is therefore "regularized_admm TV-only == split_bregman on the new differences", not "== HEAD output".
4. **m-step.** Plain CG (`src/conjgrad.cpp`, `src/gpu/conjgrad.cu`), no preconditioner, warm-started from the previous x.
5. **Shrinkage.** Isotropic over the full 3x3 Jacobian (Frobenius), eps = 1e-8 inside the square root.
6. **Bregman convention.** Goldstein-Osher, `b += grad x - d`; mu is constant (not continued). The magnitude split has its own `rho`, which may grow; `b_z` is rescaled by `rho_old/rho_new`.
7. **Layout / dtype.** `VecArray<T> = std::array<Array<T>,3>`, each `Array<T>` of dims `{n1,n2,n3}`; float or double; CPU `Array`, GPU `DeviceArray`.
8. **Support mask.** `split_bregman` passes the support only to CG. The magnitude constraint uses a centred ROI block (`Roi`, same placement as `crop3d(..., SYMMETRIC)`); voxels outside are padding. Default: padding left free (absorbs out-of-view contributions); option `zero_outside` sets z = 0 there as in section 1(d).
9. **Config.** `include/config.h` (`parse_recon_params`), `include/recon_params.h` (`ReconParams`), TOML sections `[recon_params.regularizer.*]`. New keys under `combined`: `beta` (alias of `lambda_l2`), `rho`, `rho_growth`, `rho_max`, `Ms`, `estimate_Ms`, `zero_outside`.
10. **Tests.** Plain executables registered in `tests/CMakeLists.txt` (`test_regularized_admm`, `test_gpu_regularized_admm`).
11. **Memory.** `z` and `b_z` each add one `VecArray` (3 volumes). GPU adds a transient copy during the z-update. Not yet quantified at production size (`scripts/estimate_memory.py`).

### Implementation status
- Done (CPU and GPU): Tikhonov (`lambda_l2`/`beta`, system matrix only), separate `rho`, `project_sphere` (zero vector keeps previous z; initial fallback direction +x), `m0` warm start via `x0`, monotone rho continuation with `b_z` rescaling, `estimate_Ms` / fixed `Ms`, `zero_outside`.
- CPU only: `AdmmHistory` with the four residuals, objective terms (data misfit shifted by 0.5||y||^2, TV, Tikhonov), final `m` and `b_z`. The GPU history records only the magnitude residuals, rho, Ms, final `m`, `b_z`.
- Not done: driver-level output of the `|m|/Ms` map, histogram and `rho|b_z|` map; stopping rule still uses relative dx and ||x-z|| (not all four residuals); end-to-end test with the production forward operator; usage section in docs.
- Initialization-sensitivity test uses A = I, which is too well conditioned to separate the starts (all converge to the same point); the nonconvexity is not exercised by it.

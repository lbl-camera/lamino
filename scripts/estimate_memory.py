#!/usr/bin/env python3
"""
Memory footprint estimator for Tomocam reconstruction code

Mirrors the sizing rules of src/mbir.cpp / src/gpu/mbir.cu:
  - detector (PolarGrid) is padded by padded_dim()   (include/padding.h)
  - volume n1 is padded by padded_dim(); n2/n3 stay at the unpadded
    projection size, forced odd by odd_dim()
  - A^T A is applied through Toeplitz PSF kernels (include/toeplitz.h) on a
    next_fast_dim(2n-1) grid: 6 R2C kernels per PSF set (CPU: one per dataset,
    GPU: one unified set)
"""

DTYPE_SIZES = {
    'float32': 4, 'float': 4,
    'float64': 8, 'double': 8,
    'complex64': 8, 'complex<float>': 8,
    'complex128': 16, 'complex<double>': 16
}


def format_bytes(bytes_val):
    """Format bytes as human-readable string"""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if bytes_val < 1024.0:
            return f"{bytes_val:.2f} {unit}"
        bytes_val /= 1024.0
    return f"{bytes_val:.2f} PB"


def array_memory(nx, ny, nz, dtype='float32'):
    """Calculate memory for Array<T> object"""
    return nx * ny * nz * DTYPE_SIZES.get(dtype, 4)


def n_pad(n, factor):
    """include/padding.h: n_pad"""
    n2 = int(n * factor)
    return 2 * ((n2 - n) // 2)


def padded_dim(n, factor):
    """include/padding.h: padded_dim (always odd)"""
    N = n + n_pad(n, factor)
    if N % 2 == 0:
        N -= 1
    return N


def odd_dim(n):
    """include/padding.h: odd_dim (largest odd size <= n)"""
    return n - 1 if n % 2 == 0 else n


def next_fast_dim(n):
    """include/toeplitz.h: smallest 2^a * 3^b >= n"""
    best = 1
    while best < n:
        best <<= 1
    q3 = 3
    while q3 < best:
        x = q3
        while x < n:
            x <<= 1
        best = min(best, x)
        q3 *= 3
    return best


def polar_grid_memory(nangles, nrows, ncols, dtype='float32'):
    """Calculate memory for PolarGrid<T> object"""
    # 3 Arrays for x, y, z coordinates + weights w
    coords_mem = 4 * array_memory(nangles, nrows, ncols, dtype)
    # std::vector for angles
    angles_mem = nangles * DTYPE_SIZES.get(dtype, 4)
    return coords_mem + angles_mem


def psf_memory(recon_dims_padded, dtype='float32'):
    """One ToeplitzVectorOp: 6 R2C kernels on a next_fast_dim(2n-1) grid."""
    g = tuple(next_fast_dim(2 * n - 1) for n in recon_dims_padded)
    # R2C half-spectrum: g0 * g1 * (g2/2 + 1) complex values
    cplx = 'complex64' if DTYPE_SIZES.get(dtype, 4) == 4 else 'complex128'
    kernel = array_memory(g[0], g[1], g[2] // 2 + 1, cplx)
    return g, kernel, 6 * kernel


def estimate_reconstruction_memory(recon_dims, proj_dims, nangles, dtype='float32',
                                   algorithm='mbir', gpu=False, pad_factor=1.4142,
                                   ndatasets=1):
    """
    Estimate total memory for reconstruction

    Args:
        recon_dims: tuple (n1, n2, n3) - requested volume dimensions. Only n1 is
                    used: n2/n3 are derived from the projection size
        proj_dims: tuple (nrows, ncols) - projection dimensions per angle
        nangles: number of projection angles (total over all datasets)
        dtype: data type ('float32', 'float64')
        algorithm: 'mbir', 'conjgrad', 'nesterov', 'split_bregman', 'admm'
        gpu: whether using GPU acceleration
        pad_factor: padding factor (default DEFAULT_PAD_FACTOR = 1.4142)
        ndatasets: number of datasets (CPU MBIR holds one PSF set per dataset;
                   the GPU path stacks them into one unified grid)
    """
    print(f"\n{'='*60}")
    print(f"Memory Footprint Estimate - {algorithm.upper()}")
    print(f"{'='*60}")
    print(f"Requested volume: {recon_dims}")
    print(f"Projection dims: {proj_dims} × {nangles} angles ({ndatasets} dataset(s))")
    print(f"Data type: {dtype}")
    print(f"GPU acceleration: {gpu}")
    print(f"Padding factor: {pad_factor:.4f}")
    print(f"{'-'*60}")

    cdtype = 'complex64' if DTYPE_SIZES.get(dtype, 4) == 4 else 'complex128'

    # Volume: n1 padded, n2/n3 = unpadded projection size, forced odd
    vol_dims = (padded_dim(recon_dims[0], pad_factor),
                odd_dim(proj_dims[0]), odd_dim(proj_dims[1]))
    vol = array_memory(*vol_dims, dtype)  # one scalar component
    print(f"Volume component (padded to {vol_dims}): {format_bytes(vol)}")

    # Detector is padded (odd)
    proj_dims_padded = tuple(padded_dim(d, pad_factor) for d in proj_dims)
    proj_data = array_memory(nangles, *proj_dims_padded, dtype)
    print(f"Projection data (padded to {proj_dims_padded}): {format_bytes(proj_data)}")

    polar_grid = polar_grid_memory(nangles, *proj_dims_padded, dtype)
    print(f"PolarGrid geometry (padded): {format_bytes(polar_grid)}")

    # Toeplitz PSF kernels (persistent)
    psf_grid, psf_kernel, psf_set = psf_memory(vol_dims, dtype)
    n_psf_sets = 1 if gpu else ndatasets
    psf_total = n_psf_sets * psf_set
    print(f"Toeplitz PSF: grid {psf_grid}, 6 kernels × {n_psf_sets} set(s): "
          f"{format_bytes(psf_total)}")

    # sysmat() work space: 3 padded real inputs -> 3 complex spectra, + 3
    # accumulated outputs, and one real padded buffer for the C2R inverse
    psf_real = array_memory(*psf_grid, dtype)
    psf_cplx = array_memory(psf_grid[0], psf_grid[1], psf_grid[2] // 2 + 1, cdtype)
    sysmat_work = 6 * psf_cplx + psf_real
    print(f"sysmat() FFT work space: {format_bytes(sysmat_work)}")

    # Right-hand side A^T y, initial guess x0, solution x (3 components each)
    base_vectors = 3 * 3 * vol
    print(f"yT, x0, x (3 components each): {format_bytes(base_vectors)}")

    # Algorithm-specific memory
    algo_memory = 0
    print(f"\nAlgorithm-specific buffers:")

    if algorithm == 'mbir' or algorithm == 'conjgrad':
        # CG (src/conjgrad.cpp): r, p, Ap, z = clone(r) (identity precond),
        # each a 3-component vector
        cg = 4 * 3 * vol
        algo_memory = cg
        print(f"  CG vectors (r, p, Ap, z; 3 comp.): {format_bytes(cg)}")

    elif algorithm == 'nesterov':
        # momentum and extra gradient, 3 components each
        algo_memory = 2 * 3 * vol
        print(f"  Momentum & gradient: {format_bytes(algo_memory)}")

    elif algorithm == 'split_bregman':
        x_arrays = 3 * vol
        x_old_arrays = 3 * vol
        d_arrays = 9 * vol
        b_arrays = 9 * vol
        grad_x_arrays = 9 * vol
        sk_arrays = 3 * vol
        cg_arrays = 4 * 3 * vol
        algo_memory = (x_arrays + x_old_arrays + d_arrays + b_arrays +
                       grad_x_arrays + sk_arrays + cg_arrays)
        print(f"  Solution vectors (x, x_old): {format_bytes(x_arrays + x_old_arrays)}")
        print(f"  Auxiliary variables (d, b): {format_bytes(d_arrays + b_arrays)}")
        print(f"  Gradient arrays (grad_x): {format_bytes(grad_x_arrays)}")
        print(f"  Shrinkage factors (sk): {format_bytes(sk_arrays)}")
        print(f"  CG solver buffers: {format_bytes(cg_arrays)}")

    elif algorithm == 'admm':
        # regularized_admm (src/regularized_admm.cpp), all TV + const-mag on:
        x_old = 3 * vol                  # clone of x
        z_b = 2 * 3 * vol                # magnitude splitting z, dual b
        z_hist = 3 * vol                 # z_old (history / continuation)
        d_bt = 2 * 9 * vol               # TV splitting d, dual bt (3x3)
        d_old = 9 * vol                  # d_old copy in the dual update
        grad_x = 9 * vol                 # gradient of x
        sk = vol                         # TV shrinkage scale
        rhs = 3 * vol                    # right-hand side
        mask = vol                       # support mask
        cg = 4 * 3 * vol                 # inner CG (r, p, Ap, z)
        algo_memory = (x_old + z_b + z_hist + d_bt + d_old + grad_x + sk + rhs +
                       mask + cg)
        print(f"  x_old, rhs: {format_bytes(x_old + rhs)}")
        print(f"  Magnitude splitting (z, b, z_old): {format_bytes(z_b + z_hist)}")
        print(f"  TV splitting (d, bt, d_old, grad_x): "
              f"{format_bytes(d_bt + d_old + grad_x)}")
        print(f"  sk + mask: {format_bytes(sk + mask)}")
        print(f"  Inner CG vectors: {format_bytes(cg)}")

    # Adjoint NUFFT work space (backprojection of the padded projections):
    # FINUFFT typically allocates ~2-3x the output size
    nufft_memory = 3 * proj_data
    print(f"  NUFFT buffers: {format_bytes(nufft_memory)}")

    # GPU memory (if enabled): solver state lives on the device
    gpu_memory = 0
    if gpu:
        gpu_memory = base_vectors + algo_memory + psf_total + sysmat_work
        print(f"  GPU device memory (vectors + solver + PSF + work): "
              f"{format_bytes(gpu_memory)}")

    # Host total. With GPU, PSF/solver state is on the device, host keeps data
    if gpu:
        host_total = proj_data + polar_grid + base_vectors
    else:
        host_total = (proj_data + polar_grid + base_vectors + algo_memory +
                      psf_total + sysmat_work + nufft_memory)

    print(f"{'-'*60}")
    print(f"Host memory:   {format_bytes(host_total)}")
    if gpu:
        print(f"Device memory: {format_bytes(gpu_memory)}")
    total = host_total + gpu_memory
    print(f"Total estimated memory: {format_bytes(total)}")
    print(f"{'='*60}\n")

    return total


# Example usage
if __name__ == "__main__":
    print("Tomocam Memory Footprint Estimator\n")

    # Example 1: Small reconstruction
    print("\n### Example 1: Small Volume ###")
    estimate_reconstruction_memory(
        recon_dims=(25, 256, 256),
        proj_dims=(256, 256),
        nangles=141,
        dtype='float32',
        algorithm='mbir'
    )

    # Example 2: Medium reconstruction
    print("\n### Example 2: Medium Volume ###")
    estimate_reconstruction_memory(
        recon_dims=(25, 512, 512),
        proj_dims=(512, 512),
        nangles=141,
        dtype='float32',
        algorithm='mbir'
    )

    # Example 3: Large reconstruction with GPU
    print("\n### Example 3: Large Volume with GPU ###")
    estimate_reconstruction_memory(
        recon_dims=(25, 1024, 1024),
        proj_dims=(1024, 1024),
        nangles=141,
        dtype='float32',
        algorithm='mbir',
        gpu=True
    )

    # Example 4: Constant-magnitude + TV ADMM
    print("\n### Example 4: Regularized ADMM ###")
    estimate_reconstruction_memory(
        recon_dims=(25, 512, 512),
        proj_dims=(512, 512),
        nangles=141,
        dtype='float32',
        algorithm='admm',
        ndatasets=2
    )

    # Example 5: Double precision
    print("\n### Example 5: Double Precision ###")
    estimate_reconstruction_memory(
        recon_dims=(25, 1024, 1024),
        proj_dims=(1024, 1024),
        nangles=141,
        dtype='float64',
        algorithm='conjgrad'
    )

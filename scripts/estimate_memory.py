#!/usr/bin/env python3
"""
Memory footprint estimator for Tomocam reconstruction code
"""

def format_bytes(bytes_val):
    """Format bytes as human-readable string"""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if bytes_val < 1024.0:
            return f"{bytes_val:.2f} {unit}"
        bytes_val /= 1024.0
    return f"{bytes_val:.2f} PB"

def array_memory(nx, ny, nz, dtype='float32'):
    """Calculate memory for Array<T> object"""
    dtype_sizes = {
        'float32': 4, 'float': 4,
        'float64': 8, 'double': 8,
        'complex64': 8, 'complex<float>': 8,
        'complex128': 16, 'complex<double>': 16
    }
    return nx * ny * nz * dtype_sizes.get(dtype, 4)

def polar_grid_memory(nangles, nrows, ncols, dtype='float32'):
    """Calculate memory for PolarGrid<T> object"""
    # 3 Arrays for x, y, z coordinates
    coords_mem = 3 * array_memory(nangles, nrows, ncols, dtype)
    # std::vector for angles
    dtype_sizes = {'float32': 4, 'float': 4, 'float64': 8, 'double': 8}
    angles_mem = nangles * dtype_sizes.get(dtype, 4)
    return coords_mem + angles_mem

def estimate_reconstruction_memory(recon_dims, proj_dims, nangles, dtype='float32', 
                                   algorithm='mbir', gpu=False, pad_factor=1.4142):
    """
    Estimate total memory for reconstruction
    
    Args:
        recon_dims: tuple (nx, ny, nz) - reconstruction volume dimensions
        proj_dims: tuple (nrows, ncols) - projection dimensions per angle
        nangles: number of projection angles
        dtype: data type ('float32', 'float64', etc.)
        algorithm: 'mbir', 'conjgrad', 'nesterov', 'split_bregman'
        gpu: whether using GPU acceleration
        pad_factor: padding factor (default sqrt(2) = 1.4142)
    """
    print(f"\n{'='*60}")
    print(f"Memory Footprint Estimate - {algorithm.upper()}")
    print(f"{'='*60}")
    print(f"Reconstruction volume: {recon_dims}")
    print(f"Projection dims: {proj_dims} × {nangles} angles")
    print(f"Data type: {dtype}")
    print(f"GPU acceleration: {gpu}")
    print(f"Padding factor: {pad_factor:.4f}")
    print(f"{'-'*60}")
    
    # Core data structures (unpadded)
    recon_vol_unpadded = array_memory(*recon_dims, dtype)
    print(f"Reconstruction volume (unpadded): {format_bytes(recon_vol_unpadded)}")
    
    # Calculate padded dimensions
    # Based on mbir_cg.cpp logic: size_t n_pad = 2 * (static_cast<size_t>(n * padding) / 2)
    # where padding = PAD_FACTOR - 1.0
    padding = pad_factor - 1.0
    def calc_padded_size(n):
        n_pad = 2 * int((n * padding) / 2)
        return n + n_pad
    
    recon_dims_padded = tuple(calc_padded_size(d) for d in recon_dims)
    recon_vol = array_memory(*recon_dims_padded, dtype)
    print(f"Reconstruction volume (padded to {recon_dims_padded}): {format_bytes(recon_vol)}")
    
    # Projections are also padded
    proj_data_unpadded = array_memory(nangles, *proj_dims, dtype)
    print(f"Projection data (unpadded): {format_bytes(proj_data_unpadded)}")
    
    proj_dims_padded = tuple(calc_padded_size(d) for d in proj_dims)
    proj_data = array_memory(nangles, *proj_dims_padded, dtype)
    print(f"Projection data (padded to {proj_dims_padded}): {format_bytes(proj_data)}")
    
    # PolarGrid uses padded projection dimensions
    polar_grid = polar_grid_memory(nangles, *proj_dims_padded, dtype)
    print(f"PolarGrid geometry (padded): {format_bytes(polar_grid)}")
    
    # Algorithm-specific memory
    algo_memory = 0
    print(f"\nAlgorithm-specific buffers:")
    
    if algorithm == 'mbir':
        # Gradient, residual, preconditioner, and intermediate arrays
        gradient = recon_vol
        residual = proj_data
        precond = recon_vol * 0.5  # typically sparser
        temp_arrays = 3 * recon_vol  # various intermediate computations
        algo_memory = gradient + residual + precond + temp_arrays
        print(f"  Gradient: {format_bytes(gradient)}")
        print(f"  Residual: {format_bytes(residual)}")
        print(f"  Preconditioner: {format_bytes(precond)}")
        print(f"  Temporary arrays: {format_bytes(temp_arrays)}")
        
    elif algorithm == 'conjgrad':
        # CG requires: residual, search direction, and temp vectors
        r = recon_vol
        p = recon_vol
        Ap = recon_vol
        z = recon_vol  # preconditioned residual
        algo_memory = r + p + Ap + z
        print(f"  CG vectors (r,p,Ap,z): {format_bytes(algo_memory)}")
        
    elif algorithm == 'nesterov':
        # Nesterov requires momentum and extra gradient
        momentum = recon_vol
        extra_grad = recon_vol
        algo_memory = momentum + extra_grad
        print(f"  Momentum & gradient: {format_bytes(algo_memory)}")
        
    elif algorithm == 'split_bregman':
        # Split-Bregman for vector field (3 components)
        # x, x_old: 3 arrays each
        x_arrays = 3 * recon_vol
        x_old_arrays = 3 * recon_vol
        # d, b: 3x3 arrays each (auxiliary and Bregman variables)
        d_arrays = 9 * recon_vol
        b_arrays = 9 * recon_vol
        # grad_x: 3x3 arrays (gradient of x)
        grad_x_arrays = 9 * recon_vol
        # sk: 3 arrays (shrinkage factors)
        sk_arrays = 3 * recon_vol
        # rhs, temp arrays for CG solver
        cg_arrays = 4 * 3 * recon_vol  # CG needs ~4 arrays per component
        
        algo_memory = (x_arrays + x_old_arrays + d_arrays + b_arrays + 
                      grad_x_arrays + sk_arrays + cg_arrays + proj_data)
        print(f"  Solution vectors (x, x_old): {format_bytes(x_arrays + x_old_arrays)}")
        print(f"  Auxiliary variables (d, b): {format_bytes(d_arrays + b_arrays)}")
        print(f"  Gradient arrays (grad_x): {format_bytes(grad_x_arrays)}")
        print(f"  Shrinkage factors (sk): {format_bytes(sk_arrays)}")
        print(f"  CG solver buffers: {format_bytes(cg_arrays)}")
        print(f"  Residual/projections: {format_bytes(proj_data)}")
    
    # FFT/NUFFT working memory (FINUFFT typically allocates ~2-3x output size)
    fft_memory = 3 * proj_data
    print(f"  FFT/NUFFT buffers: {format_bytes(fft_memory)}")
    
    # GPU memory (if enabled)
    gpu_memory = 0
    if gpu:
        # Approximate: device arrays for reconstruction + projections
        gpu_memory = recon_vol + proj_data + algo_memory * 0.5
        print(f"  GPU device memory: {format_bytes(gpu_memory)}")
    
    # Total
    total = recon_vol + proj_data + polar_grid + algo_memory + fft_memory + gpu_memory
    
    print(f"{'-'*60}")
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
        proj_dims=(512, 512),
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
    
    # Example 4: Double precision
    print("\n### Example 4: Double Precision ###")
    estimate_reconstruction_memory(
        recon_dims=(25, 1024, 1024),
        proj_dims=(1024, 1024),
        nangles=141,
        dtype='float64',
        algorithm='conjgrad'
    )

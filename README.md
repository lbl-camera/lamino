# Tomocam

A C++/CUDA library for reconstruction of magnetic field in materials exhibiting magetic circular dichroism (XMCD). The code is optimized for samples with thin form-factor.

[![Docker image](https://github.com/lbl-camera/lamino/actions/workflows/docker.yml/badge.svg?branch=master&event=status)](https://github.com/lbl-camera/lamino/actions/workflows/docker.yml)


## Overview

Tomocam is a high-performance library developed at Lawrence Berkeley National Laboratory for advanced reconstruction of magnetic field in thin XMCD materials. It provides forward and backward projection operators, iterative reconstruction algorithms, and Model-Based Iterative Reconstion.

The library is optimized for performance using:
- OpenMP parallelization
- Intel TBB (Threading Building Blocks)
- FFTW for fast Fourier transforms
- FINUFFT for non-uniform FFTs
- (Optional) GPU acceleration via CUDA
- (Optional) cuFINUFFT for GPU-accelerated non-uniform FFTs

## Features

- **Forward/Backward Projection**: NUFFT based projection operators
- **Iterative Reconstruction**: Conjugate gradient, Split-Bregman, and Nesterov accelerated gradient methods
- **TIFF I/O**: Read and write reconstruction data
- **Paraview**: Plot magnetization vectors in paraview

## Requirements

### Build Dependencies

- CMake ≥ 3.20
- C++ compiler with C++20 support (llvm recommended)
- OpenMP
- Intel TBB
- FFTW3 (single and double precision)
- libtiff
- FINUFFT
- CUDA (Optional)
    - cuFINUFFT
    - cufft
    - thrust


## Building

### Configure and Build

```bash
# Linux (CPU only)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build 
```

### Build with CUDA Support

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_CUDA:BOOL=ON
cmake --build build
```

## Usage

### Running Reconstruction

The reconstruction tool uses a TOML configuration file to specify input data and parameters:

```bash
./build/recon <config.toml>
```

### TOML Configuration File

Create a TOML file with the following structure for vector magnetic field reconstruction:

```toml
# Multiple input datasets with different gamma angles for vector reconstruction
[[input]]
filename = "/path/to/gamma0_stack.tiff"
angles = "/path/to/gamma0_angles.txt"
gamma = 0

[[input]]
filename = "/path/to/gamma45_stack.tiff"
angles = "/path/to/gamma45_angles.txt"
gamma = 45

[output]
filename = "output.tiff"                # Output filename for reconstruction
formats = ["tiff", "vti"]              # Available: "tiff", "vti"

[recon_params]
max_outer_iters = 50                    # Maximum outer iterations
tol = 1e-5                              # Convergence tolerance
xtol = 1e-5                             # X-tolerance for convergence
recon_dims = [51, 511, 511]            # Reconstruction dimensions [thickness, height, width]

[recon_params.regularizer]
method = "split_bregman"                # Regularizer: "split_bregman" or "qGGMRF"

# Parameters for split_bregman method
[recon_params.regularizer.split_bregman]
lambda = 0.5                            # Regularization parameter
mu = 10.0                               # Penalty parameter

```

**Notes:**
- Use multiple `[[input]]` sections to specify datasets at different gamma angles for vector field reconstruction
- Angles can be in degrees or radians (automatically converted)
- The angles file should contain one angle per line
- Run `./build/recon` without arguments to generate a template configuration file (`config.toml`)

## Documentation

📚 **Complete documentation is available on ReadTheDocs:**

- **Latest Documentation**: https://camera-lamino.readthedocs.io/
- **Getting Started Guide**: https://camera-lamino.readthedocs.io/en/latest/getting_started.html
- **Installation Instructions**: https://camera-lamino.readthedocs.io/en/latest/installation.html
- **Usage Guide**: https://camera-lamino.readthedocs.io/en/latest/usage.html
- **Configuration Reference**: https://camera-lamino.readthedocs.io/en/latest/configuration.html
- **API Documentation**: https://camera-lamino.readthedocs.io/en/latest/api/index.html
- **Examples**: https://camera-lamino.readthedocs.io/en/latest/examples.html

The documentation is automatically built from the `docs/` directory and updates with every commit.

## License

Copyright (c) 2018, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.

If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.

This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit other to do so.

## Contributing

Contributions are welcome! Please see our [Contributing Guide](https://camera-lamino.readthedocs.io/en/latest/contributing.html) for details on:

- Development setup
- Code style guidelines
- Testing requirements
- Pull request process

## Contact

For questions or issues, please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.

# Tomocam

A C++ library for tomographic reconstruction of magnetic field in thin material exhibiting magetic circular dichroism (MCD).

## Overview

Tomocam is a high-performance library developed at Lawrence Berkeley National Laboratory for advanced reconstruction of magnetic field in thin MCD materials. It provides forward and backward projection operators, iterative reconstruction algorithms, and Model-Based Iterative Reconstion.

The library is optimized for performance using:
- OpenMP parallelization
- Intel TBB (Threading Building Blocks)
- FFTW for fast Fourier transforms
- FINUFFT for non-uniform FFTs
- Optional GPU acceleration via CUDA

## Features

- **Forward/Backward Projection**: Efficient projection operators on polar grids
- **Iterative Reconstruction**: Conjugate gradient and Nesterov accelerated gradient methods
- **MBIR**: Model-Based Iterative Reconstruction with QGGMRF penalty
- **TIFF I/O**: Read and write reconstruction data

## Requirements

### Build Dependencies

- CMake ≥ 3.20
- C++ compiler with C++20 support (llvm recommended)
- OpenMP
- Intel TBB
- FFTW3 (single and double precision)
- libtiff
- FINUFFT
- CUDA for GPU acceleration

### Supported Platforms

- Linux (tested on Arch Linux)
- macOS (15.5+)

## Building

### Configure and Build

```bash
# Linux
cmake --preset arch
cmake --build --preset arch

# macOS
cmake --preset macos
cmake --build --preset macos
```

### Build with Tests

```bash
# Linux with tests
cmake --preset arch -DENABLE_TESTS=ON
cmake --build --preset arch
ctest --preset arch

# macOS with tests
cmake --preset macos -DENABLE_TESTS=ON
cmake --build --preset macos
```

## Usage

### Basic Example

```cpp
#include "tomocam.h"

// Create polar grid geometry
tomocam::PolarGrid<float> grid(/* parameters */);

// Forward projection
auto projections = tomocam::forward(volume, grid);

// Backward projection
auto reconstructed = tomocam::backward(projections, grid, volume_dims);

// MBIR reconstruction
auto mbir_result = tomocam::MBIR(
    projections,
    angles,
    recon_dims,
    max_iterations,
    sigma,
    p,
    lambda,
    gamma
);
```

## API Documentation

## License

Copyright (c) 2018, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.

If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.

This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit other to do so.

## Contributing

## Contact

For questions or issues, please contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.

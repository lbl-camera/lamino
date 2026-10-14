## [1.5.0] - 2026-10-02

### 🚜 Refactor

- Enhance GPU memory and device abstraction layer
- Update core data types and array utilities
- Update polar grid CPU and GPU implementations

### 💼 Other

- Resolve merge conflicts and complete GPU integration
- Update CMake configuration and finufft module detection
- Stage in-progress GPU changes before merging master
- Move VecArray class and related functions to gpu/vec_array.h
- Use parse_recon_params() instead of ReconParams() for parameter parsing
- Replace nested DeviceArray with VecArray in TV optimization functions
- Correct calculation of non-uniform points for nufft sampling
- In-progress changes on gpu branch before alignment work
- Alignment-parameter parity (Phase A of Toeplitz GPU port)
- Port the Toeplitz-trick gradient (Phase B)
## [frozen-v1] - 2026-05-09

### 🚀 Features

- Add transpose function and improve header guards
- Enhance test suite with new tests and improvements

### 🚜 Refactor

- Rename mbir_cgls to mbir_sb and remove cgls.cpp
- Improve ramp preconditioner implementation
- Improve optimization algorithms and MBIR implementations
- Misc improvements to polar grid, padding, and build config

### 💼 Other

- Reorder parameters in residual() to match usage (yTy before gamma)
- Rename recon_params.h to config.h
- Shared_ptr for PolarGrid so that it remains in scope for the Optimizer's lifetime.
- Swap dr and dh assignments and update x coordinate calculation
- Use clone() instead of std::move() to copy x[i] to x_old[i]
- Correctly print if method is unconstrained.
- Fixed print message
- Use int to accpet negative values
## [1.0] - 2026-02-12

### 💼 Other

- Calling constructors of x, y, z

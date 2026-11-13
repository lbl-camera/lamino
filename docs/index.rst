Tomocam-lamino Documentation
=====================

Tomocam-lamino is a high-performance C++ library for tomographic reconstruction of magnetic fields in thin materials exhibiting magnetic circular dichroism (MCD). Developed at Lawrence Berkeley National Laboratory, it provides advanced reconstruction algorithms optimized for parallel computing.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   getting_started
   installation
   usage
   configuration
   examples
   contributing
   license

Overview
--------

Tomocam-lamino provides:

* **Forward/Backward Projection**: Efficient projection operators on polar grids
* **Iterative Reconstruction**: Conjugate gradient, Split-Bregman, constant-magnitude and combined constant-magnitude/TV solvers
* **TIFF I/O**: Read and write reconstruction data
* **GPU Acceleration**: Optional CUDA support for enhanced performance

Features
--------

Performance Optimizations
~~~~~~~~~~~~~~~~~~~~~~~~~

* OpenMP parallelization for multi-core CPUs
* Intel TBB (Threading Building Blocks) for task-based parallelism
* FFTW for fast Fourier transforms
* FINUFFT for non-uniform FFTs
* Optional GPU acceleration via CUDA

Supported Platforms
~~~~~~~~~~~~~~~~~~~

* Linux (tested on Arch Linux)

Quick Start
-----------

Build the project:

.. code-block:: bash

   cmake --preset linux
   cmake --build --preset linux

Run reconstruction:

.. code-block:: bash

   ./build/recon config.toml

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

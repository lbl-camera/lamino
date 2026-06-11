API Documentation
=================

Complete C++ API reference for the Tomocam library.

.. toctree::
   :maxdepth: 2
   :caption: API Modules:

   core
   reconstruction
   io
   utilities

Core Components
---------------

The core components provide fundamental data structures and operations:

* **Array**: Multi-dimensional array container
* **PolarGrid**: Polar coordinate grid geometry
* **Projection Operators**: Forward and backward projection

Reconstruction Algorithms
-------------------------

Iterative reconstruction methods:

* **MBIR**: Model-Based Iterative Reconstruction
* **Conjugate Gradient**: CG solver
* **Split Bregman**: TV regularization
* **qGGMRF**: Generalized Gaussian MRF regularization

I/O Operations
--------------

File input/output:

* **TIFF I/O**: Read and write TIFF files
* **VTI Writer**: Export VTK ImageData format
* **Configuration**: TOML configuration parsing

Utilities
---------

Helper functions and utilities:

* **FFT**: Fast Fourier Transform wrappers
* **NUFFT**: Non-uniform FFT operations
* **Timer**: Performance profiling
* **Logger**: Logging utilities

Class Reference
---------------

.. doxygenindex::
   :project: Tomocam

Main Classes
~~~~~~~~~~~~

.. doxygenclass:: tomocam::Array
   :project: Tomocam
   :members:
   :protected-members:
   :undoc-members:

.. doxygenclass:: tomocam::PolarGrid
   :project: Tomocam
   :members:
   :protected-members:
   :undoc-members:

Function Reference
------------------

Projection Functions
~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: tomocam::forward
   :project: Tomocam

.. doxygenfunction:: tomocam::backward
   :project: Tomocam

Reconstruction Functions
~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: tomocam::MBIR
   :project: Tomocam

.. doxygenfunction:: tomocam::conjugate_gradient
   :project: Tomocam

Type Definitions
----------------

.. doxygentypedef:: tomocam::real_t
   :project: Tomocam

.. doxygentypedef:: tomocam::complex_t
   :project: Tomocam

Usage Example
-------------

Basic API usage:

.. code-block:: cpp

   #include "tomocam.h"

   int main() {
       using namespace tomocam;

       // Create polar grid
       PolarGrid<float> grid(
           /* n_radial */ 256,
           /* n_angular */ 512,
           /* radius */ 1.0f
       );

       // Load projection data
       auto projections = read_tiff<float>("projections.tiff");
       auto angles = read_angles("angles.txt");

       // Perform MBIR reconstruction
       Array3D<float> recon_dims = {32, 256, 256};
       auto result = MBIR(
           projections,
           angles,
           recon_dims,
           /* max_iter */ 50,
           /* sigma */ 1000.0f,
           /* p */ 1.2f,
           /* lambda */ 0.1f,
           /* gamma */ 0.0f
       );

       // Save result
       write_tiff("reconstruction.tiff", result);

       return 0;
   }

Build Configuration
-------------------

To use Tomocam in your project:

CMake
~~~~~

.. code-block:: cmake

   find_package(Tomocam REQUIRED)
   target_link_libraries(your_target PRIVATE tomocam::tomocam)

Compiler Flags
~~~~~~~~~~~~~~

.. code-block:: bash

   g++ -std=c++20 -I/path/to/tomocam/include \
       -L/path/to/tomocam/lib -ltomocam \
       -fopenmp -ltbb -lfftw3f -lfftw3 \
       your_code.cpp -o your_program

See Also
--------

* :doc:`../usage` - Usage guide and examples
* :doc:`../configuration` - Configuration file reference
* :doc:`../examples` - Practical examples

Installation
============

This guide covers the installation of Tomocam-lamino and its dependencies.

Requirements
------------

Build Dependencies
~~~~~~~~~~~~~~~~~~

* **CMake** ≥ 3.20
* **C++ compiler** with C++20 support (LLVM/Clang recommended)
* **OpenMP** - Parallel programming support
* **Intel TBB** - Threading Building Blocks
* **FFTW3** - Fast Fourier Transform library (single and double precision)
* **libtiff** - TIFF image format support
* **FINUFFT** - Non-uniform FFT library
* **CUDA** (optional) - GPU acceleration

Supported Platforms
~~~~~~~~~~~~~~~~~~~

* Linux (tested on Arch Linux, Ubuntu 20.04+)
* macOS 15.5 and later

Installing Dependencies
-----------------------

Linux (Arch Linux)
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   sudo pacman -S cmake clang openmp tbb fftw libtiff finufft

Linux (Ubuntu/Debian)
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   sudo apt-get update
   sudo apt-get install cmake build-essential libomp-dev libtbb-dev \
                        libfftw3-dev libtiff-dev

   # FINUFFT must be built from source
   git clone https://github.com/flatironinstitute/finufft.git
   cd finufft
   make lib
   sudo make install

macOS
~~~~~

.. code-block:: bash

   brew install cmake llvm libomp tbb fftw libtiff finufft

Building Tomocam-lamino
----------------

Clone the Repository
~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone <repository-url>
   cd lamino

Configure and Build (Linux)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   cmake --preset arch
   cmake --build --preset arch

Configure and Build (macOS)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   cmake --preset macos
   cmake --build --preset macos

Custom Build Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~

For custom builds, you can specify options directly:

.. code-block:: bash

   mkdir build && cd build
   cmake .. -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_CXX_COMPILER=clang++ \
            -DENABLE_TESTS=ON
   cmake --build . -j$(nproc)

Building with Tests
-------------------

To build and run tests:

.. code-block:: bash

   # Linux
   cmake --preset arch -DENABLE_TESTS=ON
   cmake --build --preset arch
   ctest --preset arch

   # macOS
   cmake --preset macos -DENABLE_TESTS=ON
   cmake --build --preset macos

Installation
------------

Install to system directories:

.. code-block:: bash

   sudo cmake --install build

Or specify a custom prefix:

.. code-block:: bash

   cmake --install build --prefix /opt/tomocam

This will install:

* Executables: ``recon``, ``forward`` to ``<prefix>/bin``
* Shared library: ``libtomocam.so`` to ``<prefix>/lib``

GPU Support (Optional)
----------------------

To enable CUDA GPU acceleration:

1. Install CUDA Toolkit (11.0 or later)
2. Configure with CUDA enabled:

.. code-block:: bash

   cmake .. -DENABLE_CUDA=ON
   cmake --build .

Verifying Installation
----------------------

Test the installation:

.. code-block:: bash

   ./build/recon

This should display usage information if the build was successful.

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**FINUFFT not found**

Ensure FINUFFT is installed and CMake can find it. You may need to set:

.. code-block:: bash

   export CMAKE_PREFIX_PATH=/path/to/finufft:$CMAKE_PREFIX_PATH

**OpenMP not available**

Install OpenMP support for your compiler:

.. code-block:: bash

   # For Clang on macOS
   brew install libomp

**C++20 features not supported**

Update your compiler to a recent version:

* GCC ≥ 11
* Clang ≥ 14
* MSVC ≥ 19.29

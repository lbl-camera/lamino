Getting Started
===============

This guide will help you get started with Tomocam for magnetic field reconstruction.

What is Tomocam?
----------------

Tomocam is a specialized library for reconstructing 3D magnetic field distributions from tomographic measurements of materials exhibiting magnetic circular dichroism (MCD). It implements state-of-the-art iterative reconstruction algorithms optimized for high-performance computing.

Prerequisites
-------------

Before installing Tomocam, ensure you have:

* C++ compiler with C++20 support (LLVM/Clang recommended)
* CMake ≥ 3.20
* Basic understanding of tomographic reconstruction
* Familiarity with TOML configuration files

Key Concepts
------------

Magnetic Circular Dichroism (MCD)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MCD is a spectroscopic technique that measures the difference in absorption of left and right circularly polarized light by a material in a magnetic field. This property is used to map magnetic field distributions in thin materials.

Polar Grid Geometry
~~~~~~~~~~~~~~~~~~~

Tomocam uses a polar grid system for efficient projection operations. The geometry is optimized for the specific acquisition patterns used in MCD imaging.

Reconstruction Methods
~~~~~~~~~~~~~~~~~~~~~~

The library implements several reconstruction approaches:

* **Filtered Back Projection**: Fast analytical reconstruction
* **Conjugate Gradient (CG)**: Iterative least-squares reconstruction
* **Model-Based Iterative Reconstruction (MBIR)**: Advanced reconstruction with regularization
* **Split Bregman**: Efficient algorithm for TV regularization

Vector Field Reconstruction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For complete magnetic field reconstruction, multiple datasets acquired at different gamma angles are combined to reconstruct all three components of the magnetic field vector.

First Steps
-----------

1. :doc:`Install Tomocam <installation>` on your system
2. Prepare your input data (TIFF stacks and angle files)
3. Create a :doc:`configuration file <configuration>`
4. Run the reconstruction
5. Analyze the output

Next Steps
----------

* Learn about :doc:`configuration options <configuration>`
* Explore :doc:`usage examples <examples>`
* Read the :doc:`API documentation <api/index>`

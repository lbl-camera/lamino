Usage Guide
===========

This guide covers the basic usage of Tomocam for reconstruction tasks.

Command Line Interface
----------------------

Reconstruction Tool
~~~~~~~~~~~~~~~~~~~

The main reconstruction tool reads a TOML configuration file:

.. code-block:: bash

   ./build/recon <config.toml>

Generate a template configuration:

.. code-block:: bash

   ./build/recon

This creates a ``config.toml`` template in the current directory.

Forward Projection Tool
~~~~~~~~~~~~~~~~~~~~~~~

Generate synthetic projection data from a volume:

.. code-block:: bash

   ./build/forward <input_volume.tiff> <output_projections.tiff> <angles.txt>

Workflow Overview
-----------------

1. **Prepare Input Data**
   
   * TIFF stack(s) of projection images
   * Text file with projection angles (one per line)
   * Optional: Multiple stacks at different gamma angles for vector reconstruction

2. **Create Configuration File**

   Define reconstruction parameters in TOML format

3. **Run Reconstruction**

   Execute the recon tool with your configuration

4. **Analyze Results**

   Output TIFF file contains the reconstructed volume

Input Data Format
-----------------

TIFF Stacks
~~~~~~~~~~~

* Multi-page TIFF file
* Each page is one projection image
* Supported formats: 8-bit, 16-bit, 32-bit float
* All images must have the same dimensions

Angle Files
~~~~~~~~~~~

Plain text file with one angle per line:

.. code-block:: text

   0.0
   1.5
   3.0
   4.5
   ...

Angles can be in degrees or radians (automatically detected).

Reconstruction Types
--------------------

Scalar Reconstruction
~~~~~~~~~~~~~~~~~~~~~

Single input dataset for one component:

.. code-block:: toml

   [[input]]
   filename = "projections.tiff"
   angles = "angles.txt"
   gamma = 0

   [output]
   filename = "reconstruction.tiff"

Vector Field Reconstruction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Multiple datasets at different gamma angles:

.. code-block:: toml

   [[input]]
   filename = "gamma0.tiff"
   angles = "angles0.txt"
   gamma = 0

   [[input]]
   filename = "gamma45.tiff"
   angles = "angles45.txt"
   gamma = 45

   [[input]]
   filename = "gamma90.tiff"
   angles = "angles90.txt"
   gamma = 90

   [output]
   filename = "vector_recon.tiff"

Basic Example
-------------

Complete workflow example:

.. code-block:: bash

   # 1. Generate template configuration
   ./build/recon

   # 2. Edit config.toml with your parameters
   vim config.toml

   # 3. Run reconstruction
   ./build/recon config.toml

   # 4. View results
   # The output TIFF can be opened with ImageJ, FIJI, or similar tools

Advanced Usage
--------------

Parallel Execution
~~~~~~~~~~~~~~~~~~

Control OpenMP threads:

.. code-block:: bash

   export OMP_NUM_THREADS=16
   ./build/recon config.toml

Memory Management
~~~~~~~~~~~~~~~~~

For large reconstructions, monitor memory usage:

.. code-block:: bash

   # Estimate memory requirements
   python estimate_memory.py --dims 512 512 128

Batch Processing
~~~~~~~~~~~~~~~~

Process multiple datasets:

.. code-block:: bash

   #!/bin/bash
   for config in configs/*.toml; do
       echo "Processing $config"
       ./build/recon "$config"
   done

Output Formats
--------------

TIFF Output
~~~~~~~~~~~

* Default output format
* Multi-page TIFF for 3D volumes
* 32-bit float precision
* Compatible with ImageJ, FIJI, Python (tifffile, scikit-image)

VTI Output
~~~~~~~~~~

* VTK ImageData format
* For visualization in ParaView, VisIt
* Includes spatial metadata

Enable in configuration:

.. code-block:: toml

   [output]
   formats = ["tiff", "vti"]

Visualization
-------------

ImageJ/FIJI
~~~~~~~~~~~

.. code-block:: bash

   fiji reconstruction.tiff

Python
~~~~~~

.. code-block:: python

   import tifffile
   import matplotlib.pyplot as plt

   # Load reconstruction
   volume = tifffile.imread('reconstruction.tiff')

   # View middle slice
   plt.imshow(volume[volume.shape[0]//2], cmap='gray')
   plt.colorbar()
   plt.show()

ParaView (for VTI files)
~~~~~~~~~~~~~~~~~~~~~~~~

1. Open ParaView
2. File → Open → Select .vti file
3. Apply
4. Choose visualization (Volume Rendering, Slices, etc.)

Performance Tips
----------------

* Use appropriate ``recon_dims`` - larger dimensions require more memory/time
* Adjust ``max_outer_iters`` based on convergence monitoring
* Use ``tol`` and ``xtol`` to control when reconstruction stops
* Enable GPU acceleration for large problems (if available)
* Use multiple CPU cores (automatically done with OpenMP)

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**Out of Memory**

Reduce reconstruction dimensions or process on a machine with more RAM.

**Slow Convergence**

* Increase ``max_outer_iters``
* Adjust regularization parameters
* Check input data quality

**Poor Reconstruction Quality**

* Verify angle file matches projection stack
* Check for alignment issues in input data
* Tune regularization parameters
* Ensure sufficient angular coverage

Next Steps
----------

* Learn about :doc:`configuration options <configuration>` in detail
* Explore :doc:`examples` for specific use cases
* Check the :doc:`API documentation <api/index>` for library usage

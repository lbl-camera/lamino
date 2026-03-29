Configuration Reference
=======================

Complete reference for TOML configuration files used by Tomocam-lamino

Configuration File Structure
-----------------------------

The configuration file uses TOML format with three main sections:

* ``[[input]]`` - Input data specifications, repeat for each orientation
* ``[output]`` - Output file settings
* ``[recon_params]`` - Reconstruction parameters

Input Section
-------------

Each ``[[input]]`` block specifies one dataset:

.. code-block:: toml

   [[input]]
   filename = "/path/to/projections.tiff"
   angles = "/path/to/angles.txt"
   gamma = 45 

Parameters
~~~~~~~~~~

**filename** (string, required)
   Path to TIFF stack containing projection images.

**angles** (string, required)
   Path to text file with projection angles (one per line).
   Angles can be in degrees or radians.

**gamma** (number, required)
   Orientation angle in degrees for this dataset.

Multiple Inputs
~~~~~~~~~~~~~~~

For vector reconstruction, specify multiple input blocks:

.. code-block:: toml

   [[input]]
   filename = "gamma0_stack.tiff"
   angles = "gamma0_angles.txt"
   gamma = 0

   [[input]]
   filename = "gamma45_stack.tiff"
   angles = "gamma45_angles.txt"
   gamma = 45

   [[input]]
   filename = "gamma90_stack.tiff"
   angles = "gamma90_angles.txt"
   gamma = 90

Output Section
--------------

Specifies output file settings:

.. code-block:: toml

   [output]
   filename = "reconstruction.tiff"
   formats = ["tiff", "vti"]

Parameters
~~~~~~~~~~

**filename** (string, required)
   Base filename for output. Extension is added based on format.

**formats** (array of strings, optional)
   List of output formats to generate. Default: ``["tiff"]``
   
   Available formats:
   
   * ``"tiff"`` - Multi-page TIFF file (32-bit float)
   * ``"vti"`` - VTK ImageData format for ParaView/VisIt

Reconstruction Parameters
-------------------------

Main reconstruction settings:

.. code-block:: toml

   [recon_params]
   max_outer_iters = 100
   tol = 1e-5
   xtol = 1e-5
   recon_dims = [21, 511, 511]

Parameters
~~~~~~~~~~

**max_outer_iters** (integer, required)
   Maximum number of outer iterations. Range: 1-1000.
   Typical values: 20-100.

**tol** (float, required)
   Convergence tolerance for objective function.
   Reconstruction stops when relative change < tol.
   Typical values: 1e-5 to 1e-3.

**xtol** (float, required)
   Convergence tolerance for solution update.
   Reconstruction stops when relative change in x < xtol.
   Typical values: 1e-5 to 1e-3.

**recon_dims** (array of 3 integers, required)
   Reconstruction volume dimensions ``[thickness, height, width]``.
   Tomocam-lamino will truncate even dimensions by 1.
   
   * First dimension (small compared to x- and y-axis): Sample thickness (z-axis)
   * Second dimension: Height (y-axis)
   * Third dimension: Width (x-axis)

Regularization
--------------

Tomocam-lamino supports two regularization methods:

Split Bregman Method
~~~~~~~~~~~~~~~~~~~~

Efficient total variation regularization:

.. code-block:: toml

   [recon_params.regularizer]
   method = "split_bregman"

   [recon_params.regularizer.split_bregman]
   inner_iters = 3
   lambda = 1.0
   mu = 10.0

Parameters:

**lambda** (float)
   Regularization strength. Higher values = more smoothing.
   Range: 0.001 - 10.0. Typical: 0.01 - 0.5.

**mu** (float)
   Penalty parameter for split Bregman algorithm.
   Range: 1.0 - 100.0. Typical: 5.0 - 20.0.

**inner_iters** (int)
  Split-Bregman runs a CG for data-fidelity part

QGGMRF Regularization
~~~~~~~~~~~~~~~~~~~~~

q-Generalized Gaussian Markov Random Field:

.. code-block:: toml

   [recon_params.regularizer]
   method = "qGGMRF"

   [recon_params.regularizer.qGGMRF]
   sigma = 1000.0
   p = 1.2

Parameters:

**sigma** (float)
   Scale parameter. Controls regularization strength.
   Range: 10.0 - 10000.0. Typical: 100.0 - 2000.0.

**p** (float)
   Shape parameter. Controls edge preservation.
   Range: 1.0 - 2.0. 
   
   * p = 2.0: smoothest
   * p = 1.0: preserves edges

Complete Example
----------------

Full configuration file example:

.. code-block:: toml

   # Tomocam-lamino Reconstruction Configuration

   # Input datasets - multiple for vector reconstruction
   [[input]]
   filename = "/data/mcd/gamma0_projections.tiff"
   angles = "/data/mcd/gamma0_angles.txt"
   gamma = 0

   [[input]]
   filename = "/data/mcd/gamma45_projections.tiff"
   angles = "/data/mcd/gamma45_angles.txt"
   gamma = 45

   [[input]]
   filename = "/data/mcd/gamma90_projections.tiff"
   angles = "/data/mcd/gamma90_angles.txt"
   gamma = 90

   # Output configuration
   [output]
   filename = "mcd_vector_recon.tiff"
   formats = ["tiff", "vti"]

   # Reconstruction parameters
   [recon_params]
   max_outer_iters = 100
   tol = 1e-5
   xtol = 1e-5
   recon_dims = [21, 511, 511]

   # Regularization - Split Bregman
   [recon_params.regularizer]
   method = "split_bregman"

   [recon_params.regularizer.split_bregman]
   lambda = 1.0
   mu = 10.0
   inner_iters = 3

Parameter Selection Guide
--------------------------

Choosing Reconstruction Dimensions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Consider:

* **Memory**: Larger dimensions need more RAM (~4 bytes per voxel)
* **Resolution**: Match your detector pixel size
* **Coverage**: Ensure all sample features are captured

Typical values:

* Small samples:  11 x 255 x 255
* Large samples: 21 x 511 x 511

Choosing Regularization Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**For Split Bregman:**

1. Start with ``lambda = 1.0``, ``mu = 10.0``
2. Increase lambda if reconstruction is noisy
3. Decrease lambda if reconstruction is over-smoothed
4. Adjust mu if convergence is slow (increase) or unstable (decrease)

**For qGGMRF:**

1. Start with ``sigma = 1000.0``, ``p = 1.2``
2. Increase sigma for more regularization
3. Decrease p to preserve sharper edges (minimum 1.0)
4. Increase p for smoother reconstruction (maximum 2.0)

Convergence Parameters
~~~~~~~~~~~~~~~~~~~~~~

* Set ``max_outer_iters`` high enough (50-100) for convergence
* Use ``tol`` and ``xtol`` around 1e-5 for good quality
* Monitor convergence in output logs
* Reduce tolerances if reconstruction quality is insufficient

Best Practices
--------------

1. **Start Simple**: Begin with default parameters
2. **Iterate**: Adjust one parameter at a time
3. **Validate**: Check reconstruction quality visually
4. **Document**: Keep notes on successful parameter combinations
5. **Batch**: Use same settings for similar samples

Troubleshooting
---------------

**Reconstruction too noisy**
   Increase regularization (lambda or sigma)

**Reconstruction too smooth/blurry**
   Decrease regularization

**Slow convergence**
   * Increase max_outer_iters
   * Adjust mu (split Bregman)
   * Check input data quality

**Out of memory**
   Reduce recon_dims

**Poor reconstruction quality**
   * Verify input data and angles
   * Adjust regularization method
   * Ensure sufficient angular coverage

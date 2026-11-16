Examples
========

Practical examples for common reconstruction tasks.

Example 1: Vector Field Reconstruction
---------------------------------------

Reconstruct three magnetic field components using multiple gamma orientations.

Configuration File
~~~~~~~~~~~~~~~~~~

.. code-block:: toml

   # vector_recon.toml
   [[input]]
   filename = "data/magnetic_sample_gamma0.tiff"
   angles = "data/angles.txt"
   gamma = 0

   [[input]]
   filename = "data/magnetic_sample_gamma45.tiff"
   angles = "data/angles.txt"
   gamma = 45

   [[input]]
   filename = "data/magnetic_sample_gamma90.tiff"
   angles = "data/angles.txt"
   gamma = 90

   [output]
   filename = "output/magnetic_vector_field.tiff"
   formats = ["tiff", "vti"]

   [recon_params]
   max_iters = 100
   tol = 1e-5
   xtol = 1e-5
   recon_dims = [21, 511, 511]

   [recon_params.regularizer]
   method = "split_bregman"

   [recon_params.regularizer.split_bregman]
   lambda = 1.0
   mu = 20.0

Run Reconstruction
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   export OMP_NUM_THREADS=16  # Use 16 CPU cores
   ./build/recon vector_recon.toml

Example 2: Constant Magnitude with TV
--------------------------------------

Vector reconstruction of a sample with uniform magnetization magnitude, started
from a CG data-fit solution.

Configuration File
~~~~~~~~~~~~~~~~~~

.. code-block:: toml

   # constmag_tv.toml
   [[input]]
   filename = "data/gamma0.tiff"
   angles = "data/angles.txt"
   gamma = 0

   [[input]]
   filename = "data/gamma45.tiff"
   angles = "data/angles.txt"
   gamma = 45

   [output]
   filename = "output/constmag_tv.tiff"
   formats = ["tiff", "vti"]

   [recon_params]
   max_iters = 500
   tol = 1e-5
   xtol = 1e-5
   recon_dims = [9, 319, 319]
   warm_start_iters = 50      # CG initialisation (default)

   [recon_params.regularizer]
   method = "combined"

   [recon_params.regularizer.combined]
   const_mag = true
   lambda_tv = 0.2
   mu = 10.0
   rho = 10.0
   rho_growth = 1.05
   rho_max = 200.0
   inner_iters = 5

Run Reconstruction
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Estimate memory requirements first
   python scripts/estimate_memory.py --help

   ./build/recon constmag_tv.toml

Example 3: Batch Processing
----------------------------

Process multiple samples with a script.

Batch Script
~~~~~~~~~~~~

.. code-block:: bash

   #!/bin/bash
   # batch_recon.sh

   # Create output directory
   mkdir -p batch_output

   # Sample names
   SAMPLES=("sample1" "sample2" "sample3" "sample4")

   # Process each sample
   for sample in "${SAMPLES[@]}"; do
       echo "Processing $sample..."
       
       # Create config from template
       cat > "config_${sample}.toml" << EOF
   [[input]]
   filename = "data/${sample}_projections.tiff"
   angles = "data/${sample}_angles.txt"
   gamma = 0

   [output]
   filename = "batch_output/${sample}_recon.tiff"
   formats = ["tiff"]

   [recon_params]
   max_iters = 50
   tol = 1e-5
   xtol = 1e-5
   recon_dims = [32, 256, 256]

   [recon_params.regularizer]
   method = "split_bregman"

   [recon_params.regularizer.split_bregman]
   lambda = 0.1
   mu = 10.0
   EOF

       # Run reconstruction
       ./build/recon "config_${sample}.toml"
       
       # Clean up config
       rm "config_${sample}.toml"
       
       echo "Completed $sample"
   done

   echo "All samples processed!"

Run Batch
~~~~~~~~~

.. code-block:: bash

   chmod +x batch_recon.sh
   ./batch_recon.sh

Example 4: Forward Projection (Simulation)
-------------------------------------------

Generate synthetic projection data from a known magnetization volume (mumax3
``.ovf`` files, one per component) with the ``forward`` tool, then reconstruct
them.

Forward Configuration
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: toml

   # forward.toml
   [input]
   data_path = "/path/to/sim"          # directory with the .ovf files
   component_1 = "m000000.ovf"         # vector OVF (or one file per component)
   component_2 = ""
   component_3 = ""

   [angles]                            # [begin, end) in degrees, num_projs steps
   begin = -70
   end = 71
   num_projs = 141

   [output]
   basedir = "/path/to/projections"
   dims = [451, 451]                   # detector size (rows, columns)

   [[projections]]
   gamma = 0
   filename = "gamma0.tiff"

   [[projections]]
   gamma = 45
   filename = "gamma45.tiff"

Instead of ``begin``/``end``/``num_projs``, ``[angles]`` may give
``filename = "angles.txt"``. Each ``[[projections]]`` entry also accepts
``beta`` and ``alignment`` (an ``alignlsq`` TOML that overrides the angles,
``gamma``, ``beta`` and shifts).

Generate Projections
~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   ./build/forward forward.toml

Reconstruct
~~~~~~~~~~~

Point the ``[[input]]`` entries of a ``recon`` configuration at the generated
TIFF stacks (one angles file with the same angles for each), then:

.. code-block:: bash

   ./build/recon recon.toml

Example 5: Python Integration
------------------------------

Use Tomocam output in Python analysis.

Load and Visualize
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import numpy as np
   import tifffile
   import matplotlib.pyplot as plt
   from mpl_toolkits.mplot3d import Axes3D

   # Load reconstruction
   recon = tifffile.imread('output/reconstruction0.tiff')  # component 0
   print(f"Reconstruction shape: {recon.shape}")
   print(f"Value range: [{recon.min():.3f}, {recon.max():.3f}]")

   # Visualize middle slices
   fig, axes = plt.subplots(1, 3, figsize=(15, 5))

   # XY slice (middle z)
   axes[0].imshow(recon[recon.shape[0]//2], cmap='RdBu')
   axes[0].set_title('XY Slice (middle)')
   axes[0].axis('off')

   # XZ slice (middle y)
   axes[1].imshow(recon[:, recon.shape[1]//2], cmap='RdBu')
   axes[1].set_title('XZ Slice (middle)')
   axes[1].axis('off')

   # YZ slice (middle x)
   axes[2].imshow(recon[:, :, recon.shape[2]//2], cmap='RdBu')
   axes[2].set_title('YZ Slice (middle)')
   axes[2].axis('off')

   plt.tight_layout()
   plt.savefig('slices.png', dpi=150)
   plt.show()

Compute Statistics
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import numpy as np
   import tifffile

   recon = tifffile.imread('output/reconstruction0.tiff')  # component 0

   # Basic statistics
   print(f"Mean: {np.mean(recon):.6f}")
   print(f"Std Dev: {np.std(recon):.6f}")
   print(f"Min: {np.min(recon):.6f}")
   print(f"Max: {np.max(recon):.6f}")

   # Compute gradient magnitude
   gz, gy, gx = np.gradient(recon)
   grad_mag = np.sqrt(gz**2 + gy**2 + gx**2)
   print(f"Mean gradient magnitude: {np.mean(grad_mag):.6f}")

   # Find magnetic domains (thresholding)
   threshold = 0.5 * np.max(recon)
   domains = recon > threshold
   volume_fraction = np.sum(domains) / domains.size
   print(f"Volume fraction above threshold: {volume_fraction:.2%}")

Example 6: ParaView Visualization
----------------------------------

Visualize VTI output in ParaView.

Generate VTI Output
~~~~~~~~~~~~~~~~~~~

.. code-block:: toml

   [output]
   filename = "paraview_recon"
   formats = ["vti"]

ParaView Script
~~~~~~~~~~~~~~~

.. code-block:: python

   # paraview_script.py
   # Run in ParaView's Python shell or pvpython

   from paraview.simple import *

   # Load VTI file
   reader = XMLImageDataReader(FileName=['paraview_recon.vti'])

   # Create render view
   renderView = CreateView('RenderView')
   renderView.ViewSize = [1200, 800]

   # Volume rendering
   display = Show(reader, renderView)
   display.Representation = 'Volume'

   # Configure color map
   colorMap = GetColorTransferFunction('ImageScalars')
   colorMap.ApplyPreset('Cool to Warm', True)

   # Render
   renderView.ResetCamera()
   Render()

   # Save screenshot
   SaveScreenshot('paraview_volume.png', renderView, ImageResolution=[1200, 800])

Next Steps
----------

* Adapt these examples for your specific datasets
* Experiment with different :doc:`configuration parameters <configuration>`
* Review :doc:`usage guide <usage>` for more details

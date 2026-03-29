Examples
========

Practical examples for common reconstruction tasks.

Example 1: Basic Scalar Reconstruction
---------------------------------------

Reconstruct a single MCD component from projection data.

Configuration File
~~~~~~~~~~~~~~~~~~

.. code-block:: toml

   # basic_recon.toml
   [[input]]
   filename = "data/sample1_projections.tiff"
   angles = "data/sample1_angles.txt"
   gamma = 0

   [output]
   filename = "output/sample1_recon.tiff"
   formats = ["tiff"]

   [recon_params]
   max_outer_iters = 50
   tol = 1e-5
   xtol = 1e-5
   recon_dims = [32, 256, 256]

   [recon_params.regularizer]
   method = "split_bregman"

   [recon_params.regularizer.split_bregman]
   lambda = 0.1
   mu = 10.0

Run Reconstruction
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   ./build/recon basic_recon.toml

Example 2: Vector Field Reconstruction
---------------------------------------

Reconstruct all three magnetic field components using multiple gamma angles.

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
   max_outer_iters = 80
   tol = 1e-6
   xtol = 1e-6
   recon_dims = [64, 512, 512]

   [recon_params.regularizer]
   method = "split_bregman"

   [recon_params.regularizer.split_bregman]
   lambda = 0.05
   mu = 15.0

Run Reconstruction
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   export OMP_NUM_THREADS=16  # Use 16 CPU cores
   ./build/recon vector_recon.toml

Example 3: High-Resolution with qGGMRF
---------------------------------------

Large reconstruction with edge-preserving regularization.

Configuration File
~~~~~~~~~~~~~~~~~~

.. code-block:: toml

   # high_res_qggmrf.toml
   [[input]]
   filename = "data/highres_projections.tiff"
   angles = "data/highres_angles.txt"
   gamma = 0

   [output]
   filename = "output/highres_recon.tiff"
   formats = ["tiff", "vti"]

   [recon_params]
   max_outer_iters = 100
   tol = 1e-6
   xtol = 1e-6
   recon_dims = [128, 1024, 1024]

   [recon_params.regularizer]
   method = "qGGMRF"

   [recon_params.regularizer.qGGMRF]
   sigma = 1500.0
   p = 1.1  # Sharper edges

Run Reconstruction
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Estimate memory requirements first
   python estimate_memory.py --dims 1024 1024 128

   # Run with all available cores
   ./build/recon high_res_qggmrf.toml

Example 4: Batch Processing
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
   max_outer_iters = 50
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

Example 5: Forward Projection (Simulation)
-------------------------------------------

Generate synthetic projection data from a known volume.

Create Synthetic Volume
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # create_phantom.py
   import numpy as np
   import tifffile

   # Create a 3D phantom with magnetic features
   dims = (32, 256, 256)
   phantom = np.zeros(dims, dtype=np.float32)

   # Add sphere
   center = np.array([16, 128, 128])
   radius = 40
   for z in range(dims[0]):
       for y in range(dims[1]):
           for x in range(dims[2]):
               if np.linalg.norm([z-center[0], y-center[1], x-center[2]]) < radius:
                   phantom[z, y, x] = 1.0

   # Save phantom
   tifffile.imwrite('phantom.tiff', phantom)
   print("Phantom created: phantom.tiff")

Create Angles File
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # create_angles.py
   import numpy as np

   # Create evenly spaced angles from 0 to 180 degrees
   angles = np.linspace(0, 180, 120, endpoint=False)
   np.savetxt('angles.txt', angles, fmt='%.6f')
   print(f"Created {len(angles)} angles")

Generate Projections
~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Create phantom and angles
   python create_phantom.py
   python create_angles.py

   # Generate forward projections
   ./build/forward phantom.tiff projections.tiff angles.txt

Reconstruct
~~~~~~~~~~~

.. code-block:: bash

   # Now reconstruct from projections
   ./build/recon simulation_config.toml

Example 6: Python Integration
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
   recon = tifffile.imread('output/reconstruction.tiff')
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

   recon = tifffile.imread('output/reconstruction.tiff')

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

Example 7: ParaView Visualization
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
* Explore the :doc:`API documentation <api/index>` for custom applications
* Review :doc:`usage guide <usage>` for more details

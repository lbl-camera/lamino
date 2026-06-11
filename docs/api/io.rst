I/O API
=======

File input/output operations.

TIFF Operations
---------------

Read TIFF Files
~~~~~~~~~~~~~~~

.. code-block:: cpp

   template<typename T>
   Array3D<T> read_tiff(const std::string& filename);

Read a multi-page TIFF file into a 3D array.

Parameters:
   * **filename**: Path to TIFF file

Returns:
   3D array containing image data

Write TIFF Files
~~~~~~~~~~~~~~~~

.. code-block:: cpp

   template<typename T>
   void write_tiff(const std::string& filename, const Array3D<T>& data);

Write a 3D array to multi-page TIFF file.

Parameters:
   * **filename**: Output file path
   * **data**: 3D array to write

VTK ImageData (VTI)
-------------------

Write VTI Files
~~~~~~~~~~~~~~~

.. code-block:: cpp

   template<typename T>
   void write_vti(
       const std::string& filename,
       const Array3D<T>& data,
       const std::array<T, 3>& spacing = {1.0, 1.0, 1.0},
       const std::array<T, 3>& origin = {0.0, 0.0, 0.0}
   );

Write a 3D array to VTK ImageData XML format (.vti).

Parameters:
   * **filename**: Output file path
   * **data**: 3D array to write
   * **spacing**: Voxel spacing [x, y, z]
   * **origin**: Volume origin [x, y, z]

Configuration
-------------

TOML Configuration Parser
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: cpp

   struct ReconConfig {
       struct Input {
           std::string filename;
           std::string angles;
           float gamma;
       };
       
       struct Output {
           std::string filename;
           std::vector<std::string> formats;
       };
       
       struct ReconParams {
           size_t max_outer_iters;
           float tol;
           float xtol;
           std::array<size_t, 3> recon_dims;
           
           struct Regularizer {
               std::string method;
               // Split Bregman params
               float lambda;
               float mu;
               // qGGMRF params
               float sigma;
               float p;
           } regularizer;
       };
       
       std::vector<Input> inputs;
       Output output;
       ReconParams recon_params;
   };

   ReconConfig parse_config(const std::string& filename);

Parse a TOML configuration file.

Parameters:
   * **filename**: Path to TOML config file

Returns:
   ReconConfig structure

Angle Files
-----------

Read Angles
~~~~~~~~~~~

.. code-block:: cpp

   std::vector<float> read_angles(const std::string& filename);

Read projection angles from text file.

Parameters:
   * **filename**: Path to angles file (one angle per line)

Returns:
   Vector of angles (converted to radians if in degrees)

Write Angles
~~~~~~~~~~~~

.. code-block:: cpp

   void write_angles(
       const std::string& filename,
       const std::vector<float>& angles,
       bool degrees = false
   );

Write angles to text file.

Parameters:
   * **filename**: Output file path
   * **angles**: Vector of angles (in radians)
   * **degrees**: Write in degrees instead of radians

Example Usage
-------------

.. code-block:: cpp

   #include "tomocam.h"

   // Read input data
   auto projections = read_tiff<float>("projections.tiff");
   auto angles = read_angles("angles.txt");

   // Parse configuration
   auto config = parse_config("config.toml");

   // Perform reconstruction
   auto result = MBIR(
       projections,
       angles,
       config.recon_params.recon_dims,
       config.recon_params.max_outer_iters,
       config.recon_params.regularizer.sigma,
       config.recon_params.regularizer.p
   );

   // Write output
   write_tiff("reconstruction.tiff", result);
   write_vti("reconstruction.vti", result);

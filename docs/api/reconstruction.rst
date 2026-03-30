Reconstruction API
==================

Reconstruction algorithms and solvers.

Forward Projection
------------------

Project a 3D volume onto 2D detector plane.

.. code-block:: cpp

   template<typename T>
   Array2D<T> forward(
       const Array3D<T>& volume,
       const PolarGrid<T>& grid,
       T angle,
       T gamma = 0.0
   );

Parameters:
   * **volume**: 3D reconstruction volume
   * **grid**: Polar grid geometry
   * **angle**: Projection angle (radians)
   * **gamma**: Sample tilt angle (radians)

Returns:
   2D projection image

Backward Projection
-------------------

Backproject 2D projections into 3D volume.

.. code-block:: cpp

   template<typename T>
   Array3D<T> backward(
       const std::vector<Array2D<T>>& projections,
       const PolarGrid<T>& grid,
       const std::vector<T>& angles,
       const std::array<size_t, 3>& volume_dims,
       T gamma = 0.0
   );

Parameters:
   * **projections**: Vector of 2D projection images
   * **grid**: Polar grid geometry
   * **angles**: Projection angles (radians)
   * **volume_dims**: Output volume dimensions [z, y, x]
   * **gamma**: Sample tilt angle (radians)

Returns:
   3D reconstructed volume

MBIR Reconstruction
-------------------

Model-Based Iterative Reconstruction.

.. code-block:: cpp

   template<typename T>
   Array3D<T> MBIR(
       const std::vector<Array2D<T>>& projections,
       const std::vector<T>& angles,
       const std::array<size_t, 3>& recon_dims,
       size_t max_iterations,
       T sigma,
       T p,
       T lambda = 0.0,
       T gamma = 0.0
   );

Parameters:
   * **projections**: Input projection images
   * **angles**: Projection angles
   * **recon_dims**: Reconstruction dimensions [z, y, x]
   * **max_iterations**: Maximum iterations
   * **sigma**: qGGMRF scale parameter
   * **p**: qGGMRF shape parameter (1.0-2.0)
   * **lambda**: Regularization weight
   * **gamma**: Sample tilt angle

Returns:
   Reconstructed 3D volume

Conjugate Gradient
------------------

Conjugate gradient solver.

.. code-block:: cpp

   template<typename T>
   Array3D<T> conjugate_gradient(
       const std::function<Array3D<T>(const Array3D<T>&)>& A,
       const Array3D<T>& b,
       const Array3D<T>& x0,
       size_t max_iterations,
       T tolerance
   );

Parameters:
   * **A**: Linear operator (function)
   * **b**: Right-hand side
   * **x0**: Initial guess
   * **max_iterations**: Maximum CG iterations
   * **tolerance**: Convergence tolerance

Returns:
   Solution vector

Split Bregman
-------------

Split Bregman algorithm for TV regularization.

.. code-block:: cpp

   template<typename T>
   Array3D<T> split_bregman(
       const std::function<Array3D<T>(const Array3D<T>&)>& A,
       const std::function<Array3D<T>(const Array3D<T>&)>& AT,
       const Array3D<T>& b,
       T lambda,
       T mu,
       size_t max_outer_iters,
       size_t max_inner_iters,
       T tolerance
   );

Parameters:
   * **A**: Forward operator
   * **AT**: Adjoint operator
   * **b**: Measurement data
   * **lambda**: TV regularization parameter
   * **mu**: Penalty parameter
   * **max_outer_iters**: Outer iterations
   * **max_inner_iters**: Inner CG iterations
   * **tolerance**: Convergence tolerance

Returns:
   Reconstructed volume

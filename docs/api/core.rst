Core API
========

Core data structures and basic operations.

Array Container
---------------

Multi-dimensional array template class.

.. code-block:: cpp

   template<typename T>
   class Array {
   public:
       Array(size_t width, size_t height, size_t depth = 1);
       
       T& operator()(size_t x, size_t y, size_t z = 0);
       const T& operator()(size_t x, size_t y, size_t z = 0) const;
       
       size_t width() const;
       size_t height() const;
       size_t depth() const;
       size_t size() const;
       
       T* data();
       const T* data() const;
   };

PolarGrid
---------

Polar coordinate grid for projection geometry.

.. code-block:: cpp

   template<typename T>
   class PolarGrid {
   public:
       PolarGrid(size_t n_radial, size_t n_angular, T radius);
       
       size_t n_radial() const;
       size_t n_angular() const;
       T radius() const;
       
       // Grid coordinates
       T r(size_t i) const;  // Radial coordinate
       T theta(size_t j) const;  // Angular coordinate
   };

Data Types
----------

Type aliases for common data types:

.. code-block:: cpp

   namespace tomocam {
       using real_t = float;
       using complex_t = std::complex<float>;
       
       template<typename T>
       using Array2D = Array<T>;
       
       template<typename T>
       using Array3D = Array<T>;
   }

Utilities API
=============

Helper functions and utilities.

FFT Operations
--------------

FFT Wrapper
~~~~~~~~~~~

.. code-block:: cpp

   template<typename T>
   class FFT {
   public:
       FFT(size_t nx, size_t ny, size_t nz = 1);
       ~FFT();
       
       void forward(const T* input, std::complex<T>* output);
       void backward(const std::complex<T>* input, T* output);
       
       void forward_inplace(std::complex<T>* data);
       void backward_inplace(std::complex<T>* data);
   };

Initialize FFT plan and execute forward/backward transforms.

NUFFT Operations
----------------

Non-Uniform FFT
~~~~~~~~~~~~~~~

.. code-block:: cpp

   template<typename T>
   class NUFFT {
   public:
       NUFFT(size_t n_modes, size_t n_points);
       
       void set_points(const T* x, const T* y, const T* z = nullptr);
       
       void forward(const std::complex<T>* input, std::complex<T>* output);
       void backward(const std::complex<T>* input, std::complex<T>* output);
   };

Non-uniform FFT using FINUFFT.

Timer
-----

Performance Profiling
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: cpp

   class Timer {
   public:
       Timer();
       
       void start();
       void stop();
       
       double elapsed() const;  // seconds
       double elapsed_ms() const;  // milliseconds
       
       static Timer& global();  // Global timer instance
   };

High-resolution timer for performance measurements.

Example usage:

.. code-block:: cpp

   Timer timer;
   timer.start();
   
   // ... some computation ...
   
   timer.stop();
   std::cout << "Elapsed time: " << timer.elapsed() << " seconds\n";

Logger
------

Logging Utility
~~~~~~~~~~~~~~~

.. code-block:: cpp

   enum class LogLevel {
       DEBUG,
       INFO,
       WARNING,
       ERROR,
       CRITICAL
   };

   class Logger {
   public:
       static Logger& get();
       
       void set_level(LogLevel level);
       void set_output(std::ostream& stream);
       
       void debug(const std::string& message);
       void info(const std::string& message);
       void warning(const std::string& message);
       void error(const std::string& message);
       void critical(const std::string& message);
   };

   // Convenience macros
   #define LOG_DEBUG(msg) Logger::get().debug(msg)
   #define LOG_INFO(msg) Logger::get().info(msg)
   #define LOG_WARNING(msg) Logger::get().warning(msg)
   #define LOG_ERROR(msg) Logger::get().error(msg)

Example usage:

.. code-block:: cpp

   #include "logger.h"

   Logger::get().set_level(LogLevel::INFO);
   
   LOG_INFO("Starting reconstruction");
   LOG_DEBUG("Debug information");  // Not printed (below INFO level)
   LOG_WARNING("Convergence slow");
   LOG_ERROR("Failed to open file");

Array Operations
----------------

Basic Operations
~~~~~~~~~~~~~~~~

.. code-block:: cpp

   // Element-wise operations
   template<typename T>
   Array<T> operator+(const Array<T>& a, const Array<T>& b);
   
   template<typename T>
   Array<T> operator-(const Array<T>& a, const Array<T>& b);
   
   template<typename T>
   Array<T> operator*(const Array<T>& a, const Array<T>& b);
   
   template<typename T>
   Array<T> operator*(const Array<T>& a, T scalar);

Statistical Operations
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: cpp

   template<typename T>
   T sum(const Array<T>& a);
   
   template<typename T>
   T mean(const Array<T>& a);
   
   template<typename T>
   T std_dev(const Array<T>& a);
   
   template<typename T>
   T min(const Array<T>& a);
   
   template<typename T>
   T max(const Array<T>& a);

Gradient Operations
~~~~~~~~~~~~~~~~~~~

.. code-block:: cpp

   template<typename T>
   std::array<Array<T>, 3> gradient(const Array3D<T>& input);

Compute 3D gradient (returns [grad_z, grad_y, grad_x]).

.. code-block:: cpp

   template<typename T>
   Array3D<T> divergence(
       const Array3D<T>& dx,
       const Array3D<T>& dy,
       const Array3D<T>& dz
   );

Compute divergence of a vector field.

Padding Operations
------------------

.. code-block:: cpp

   template<typename T>
   Array<T> pad(const Array<T>& input, size_t pad_width);
   
   template<typename T>
   Array<T> unpad(const Array<T>& input, size_t pad_width);

Add or remove zero-padding around an array.

Filtering
---------

Gaussian Filter
~~~~~~~~~~~~~~~

.. code-block:: cpp

   template<typename T>
   Array<T> gaussian_filter(const Array<T>& input, T sigma);

Apply Gaussian smoothing.

Median Filter
~~~~~~~~~~~~~

.. code-block:: cpp

   template<typename T>
   Array<T> median_filter(const Array<T>& input, size_t kernel_size);

Apply median filter for noise reduction.

Memory Utilities
----------------

Memory Estimation
~~~~~~~~~~~~~~~~~

.. code-block:: cpp

   size_t estimate_memory(
       const std::array<size_t, 3>& dims,
       size_t n_projections
   );

Estimate memory requirements for reconstruction (in bytes).

.. code-block:: cpp

   bool has_sufficient_memory(size_t required_bytes);

Check if system has sufficient available memory.

Contributing
============

We welcome contributions to Tomocam! This guide will help you get started.

Getting Started
---------------

1. Fork the repository on GitHub
2. Clone your fork locally
3. Create a new branch for your feature or bugfix
4. Make your changes
5. Submit a pull request

Development Setup
-----------------

Install Development Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to the build dependencies, install:

.. code-block:: bash

   # Install clang-format for code formatting
   sudo apt-get install clang-format  # Ubuntu/Debian
   brew install clang-format          # macOS

   # Install testing tools
   pip install pytest numpy tifffile

Build with Tests
~~~~~~~~~~~~~~~~

.. code-block:: bash

   cmake --preset arch -DENABLE_TESTS=ON
   cmake --build --preset arch
   ctest --preset arch

Code Style
----------

C++ Style Guidelines
~~~~~~~~~~~~~~~~~~~~

* Follow C++20 standard
* Use modern C++ features (auto, range-based for, etc.)
* Prefer RAII and smart pointers
* Use const correctness
* Keep functions small and focused

Formatting
~~~~~~~~~~

Use clang-format with the project's .clang-format file:

.. code-block:: bash

   clang-format -i src/*.cpp include/*.h

Naming Conventions
~~~~~~~~~~~~~~~~~~

* **Classes**: PascalCase (e.g., ``PolarGrid``)
* **Functions**: snake_case (e.g., ``forward_projection``)
* **Variables**: snake_case (e.g., ``recon_dims``)
* **Constants**: UPPER_SNAKE_CASE (e.g., ``MAX_ITERATIONS``)
* **Template parameters**: PascalCase (e.g., ``typename T``)

Documentation
-------------

Code Documentation
~~~~~~~~~~~~~~~~~~

Use Doxygen-style comments:

.. code-block:: cpp

   /**
    * @brief Performs forward projection of a 3D volume
    *
    * Projects a 3D volume onto a 2D detector plane using the
    * specified polar grid geometry and projection angle.
    *
    * @param volume Input 3D volume to project
    * @param grid Polar grid geometry
    * @param angle Projection angle in radians
    * @param gamma Sample tilt angle in radians (default: 0)
    * @return 2D projection image
    */
   template<typename T>
   Array2D<T> forward(
       const Array3D<T>& volume,
       const PolarGrid<T>& grid,
       T angle,
       T gamma = 0.0
   );

Commit Messages
~~~~~~~~~~~~~~~

Write clear, descriptive commit messages:

.. code-block:: text

   Short summary (50 chars or less)

   More detailed explanation if needed. Wrap at 72 characters.
   Explain the problem this commit solves and how it solves it.

   - Use bullet points for multiple changes
   - Reference issues with #123

Testing
-------

Writing Tests
~~~~~~~~~~~~~

Add tests for new features in the ``tests/`` directory:

.. code-block:: cpp

   #include <catch2/catch.hpp>
   #include "tomocam.h"

   TEST_CASE("Forward projection", "[projection]") {
       // Arrange
       tomocam::Array3D<float> volume(256, 256, 32);
       tomocam::PolarGrid<float> grid(256, 512, 1.0f);
       
       // Act
       auto projection = tomocam::forward(volume, grid, 0.0f);
       
       // Assert
       REQUIRE(projection.width() == 256);
       REQUIRE(projection.height() == 512);
   }

Running Tests
~~~~~~~~~~~~~

.. code-block:: bash

   # Run all tests
   ctest --preset arch

   # Run specific test
   ./build/tests/test_projection

   # Run with verbose output
   ctest --preset arch --verbose

Pull Request Process
--------------------

1. **Update Documentation**: Ensure docs reflect your changes
2. **Add Tests**: Include tests for new functionality
3. **Run Tests**: Verify all tests pass
4. **Format Code**: Run clang-format on modified files
5. **Write Description**: Clearly describe what and why
6. **Reference Issues**: Link related issues in PR description

Pull Request Template
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: markdown

   ## Description
   Brief description of changes

   ## Type of Change
   - [ ] Bug fix
   - [ ] New feature
   - [ ] Performance improvement
   - [ ] Documentation update

   ## Testing
   - [ ] All existing tests pass
   - [ ] New tests added for new functionality
   - [ ] Manual testing performed

   ## Related Issues
   Fixes #123

Code Review
-----------

What We Look For
~~~~~~~~~~~~~~~~

* **Correctness**: Does it work as intended?
* **Performance**: Are there efficiency concerns?
* **Readability**: Is the code clear and well-documented?
* **Style**: Does it follow project conventions?
* **Tests**: Are changes adequately tested?

Reporting Issues
----------------

Bug Reports
~~~~~~~~~~~

Include:

* Operating system and version
* Compiler and version
* Steps to reproduce
* Expected vs actual behavior
* Relevant logs or error messages

Feature Requests
~~~~~~~~~~~~~~~~

Include:

* Use case description
* Proposed solution
* Alternative approaches considered
* Example API or usage

Community
---------

* Be respectful and constructive
* Help others learn
* Give credit where due
* Follow the code of conduct

License
-------

By contributing, you agree that your contributions will be licensed
under the same license as the project (see LICENSE file).

Questions?
----------

Contact Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov

Thank you for contributing to Tomocam!

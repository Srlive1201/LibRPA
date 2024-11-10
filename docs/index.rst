.. LibRPA documentation master file, created by
   sphinx-quickstart on Thu Sep  7 16:55:38 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

LibRPA
======

**LibRPA** is a comprehensive library designed for efficient and accurate first-principles calculations of random-phase approximation (RPA) correlation energies using numerical atomic orbitals (NAOs). It leverages the Localized Resolution of Identity (LRI) technique to achieve :math:`O(N^2)` scaling behavior, making it highly suitable for large-scale periodic systems. Implemented in C++ and Python with MPI/OpenMP parallelism, LibRPA integrates seamlessly with various density functional theory (DFT) packages through flexible file-based and API-based interfaces.

**Key Features:**

- Accurate computation of physical properties using many-body perturbation theory.
- Efficient handling of large-scale periodic systems with :math:`O(N^2)` scaling.
- Hybrid parallelism (MPI/OpenMP) to enhance scalability and performance.
- Seamless integration with various DFT packages through flexible interfaces.


**Quick Links:**

- :doc:`Installation Guide <install>`
- :doc:`Driver Usage <driver>`
- :doc:`Compile Options <user_guide/compile_options>`
- :doc:`Input Parameters <user_guide/input_parameters>`

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Get Started

   install
   driver

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: User Guide

   user_guide/compile_options
   user_guide/input_parameters
   user_guide/api_usage

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Tutorial

   tutorial/rpa/index
   tutorial/gw/index

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Developer Zone

   develop/documentation
   develop/design
   develop/dataset_format

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: API Documentation

   doxygen/librpa/html/index


**How to Cite:**

- Shi R, Lin P, Zhang M Y, et al. Subquadratic-scaling real-space random phase approximation correlation energy calculations for periodic systems with numerical atomic orbitals[J]. Physical Review B, 2024, 109(3): 035103.

..
   Indices and tables
   ==================
   
   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`

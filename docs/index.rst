.. LibRPA documentation master file, created by
   sphinx-quickstart on Thu Sep  7 16:55:38 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

LibRPA
======

**LibRPA** is a library for efficient and accurate first-principles calculations
based on many-body perturbation theory with numerical atomic orbitals (NAOs).
It currently supports random-phase approximation (RPA) correlation energy
calculations and GW quasiparticle calculations for periodic systems. By
leveraging the localized resolution of identity (LRI) technique, LibRPA
achieves low-scaling algorithms suitable for large-scale systems. Implemented
in C++ with MPI/OpenMP parallelism and equipped with Python and Fortran
interfaces, LibRPA can be integrated with density functional theory (DFT) codes
through both file-based and API-based workflows.

**Key Features**

- RPA correlation energy calculations and GW quasiparticle calculations based on NAO basis sets.
- Low-scaling algorithms enabled by the localized resolution of identity (LRI) technique.
- Efficient treatment of large-scale periodic systems.
- Hybrid MPI/OpenMP parallelism for scalability and performance.
- Flexible integration with external DFT packages through driver-based and API-based interfaces.

**Quick Links**

- :doc:`Installation Guide <install>`
- :doc:`Quick Examples <quick_examples>`
- :doc:`Compile Options <user_guide/compile_options>`
- :doc:`Runtime Parameters <user_guide/runtime_parameters>`
- :doc:`Driver Usage <user_guide/driver_usage>`
- :doc:`API Usage <user_guide/api_usage>`

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Get Started

   install
   quick_examples

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Examples

   examples/build/index
   examples/driver/index
   examples/api/index

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: User Guide

   user_guide/compile_options
   user_guide/runtime_parameters
   user_guide/driver_usage
   user_guide/api_usage

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Tutorials

   tutorial/rpa/index
   tutorial/gw/index

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Developer Zone

   develop/documentation
   develop/design
   develop/dataset_format
   develop/known_issues
   API Documentation <doxygen/librpa/html/index>


**How to Cite:**

If you use LibRPA in your work, please consider citing:

- R. Shi, P. Lin, M.-Y. Zhang, L. He, and X. Ren,
  Subquadratic-scaling real-space random phase approximation correlation energy calculations for periodic systems with numerical atomic orbitals,
  *Phys. Rev. B* **109**, 035103 (2024).

  .. code-block:: bibtex

     @article{ShiR2024,
       title = {Subquadratic-scaling real-space random phase approximation correlation energy calculations for periodic systems with numerical atomic orbitals},
       author = {Shi, Rong and Lin, Peize and Zhang, Min-Ye and He, Lixin and Ren, Xinguo},
       date = {2024-01-02},
       year = 2024,
       journal = {Phys. Rev. B},
       volume = {109},
       number = {3},
       pages = {035103},
       doi = {10.1103/PhysRevB.109.035103},
     }

- R. Shi, M.-Y. Zhang, P. Lin, L. He, and X. Ren,
  LibRPA: A software package for low-scaling first-principles calculations of random phase approximation electron correlation energy based on numerical atomic orbitals,
  *Comput. Phys. Commun.* **309**, 109496 (2025).

  .. code-block:: bibtex

     @article{ShiR25,
       title = {{{LibRPA}}: {{A}} software package for low-scaling first-principles calculations of random phase approximation electron correlation energy based on numerical atomic orbitals},
       author = {Shi, Rong and Zhang, Min-Ye and Lin, Peize and He, Lixin and Ren, Xinguo},
       date = {2025-04-01},
       year = 2025,
       journal = {Comput. Phys. Commun.},
       volume = {309},
       pages = {109496},
       doi = {10.1016/j.cpc.2024.109496},
     }

..
   Indices and tables
   ==================
   
   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`

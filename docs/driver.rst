Driver Usage
============

Data preparation
----------------

To run the LibRPA driver for RPA/GW calculation, you need to export necessary
data from atomic-basis first-principles code.
The data include

* Atomistic structure
* one-particle orbital energies, occupation numbers and wave functions
* Bloch vectors used in periodic calculation
* Triple co-efficients of resolution of identity (RI)
* Coulomb matrix under the auxiliary basis for RI

FHI-aims
~~~~~~~~

Export of data required by LibRPA from FHI-aims is supported in the latest
master branch. You can switch it on by adding to your ``control.in`` file:

.. code-block:: text

   print_input_librpa .true.

For RPA calculation (``total_energy_method rpa`` in ``control.in``), this will dump a few text files:

.. code-block:: text

   stru_out
   band_out
   KS_eigenvector_<myid>.txt
   Cs_data_<myid>.txt
   coulomb_mat_<myid>.txt

For *periodic* GW calculation (``qpe_calc gw_expt`` in ``control.in``), additional files will be created:

.. code-block:: text

   coulomb_cut_<myid>.txt

ABACUS
~~~~~~

TBA

LibRPA run
----------

The LibRPA driver can be called using the following syntax:

.. code-block:: bash

   /path/to/LibRPA/build/librpa_driver.exe <ngrids> <green-func-thres>

where ``<ngrids>`` is the number of time/frequency grids to use and
``<green-func-thres>`` is the threshold to prune the Green's function in
time-domain and real-space. To run with MPI, you can simply invoke the relavant MPI driver,
for example

.. code-block:: bash

   mpirun -np <nprocs> /path/to/LibRPA/build/librpa_driver.exe <ngrids> <green-func-thres>

where ``<nprocs>`` is the number of MPI processes.
If additionaly you want to run with multiple threads, you need to specify the
environment variable ``OMP_NUM_THREADS``. For example

.. code-block:: bash

   export OMP_NUM_THREADS=4
   mpirun -np 4 /path/to/LibRPA/build/librpa_driver.exe 16 1e-12

This will run the LibRPA RPA calculation with 16 time/frequency grids using
4 MPI processes, each with 4 OpenMP threads.


Without an input file ``librpa.in``, the above commands will all run the
low-scaling RPA calculation.
To perform other tasks such as exact-exchange or GW calculations, you
need to specify the ``task`` keyword in ``librpa.in``. Please turn to the
:doc:`manual page of input parameters <user_guide/input_parameters>` for more
information.

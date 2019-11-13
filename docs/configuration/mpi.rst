.. _configuration_mpi:

mpi
================================================================================

These parameters directly affect how UMD-LETKF scatters the model state across
the processors and how it optimizes file I/O.

.. warning::

  In an attempt to improve memory usage with the MPI calls, I over-optimized,
  resulting in overly slow performance in the MPI scatter/gather calls if a
  large domain or high number of ensemble members are used. This will be fixed
  in the future.
  
.. rubric:: Parameters
   
:ens_size:

  **type:** *integer*, ``required``  

  The number of ensemble members.


:ppn:

   **type:** *integer*,  **default:** *1*

   The number of processors per node.
   
   Optional, but may help improve I/O performance by evenly distributing across
   nodes the the PEs which perform simultaneous I/O.


.. rubric:: Example

.. code-block :: yaml

  mpi:
    ens_size: 20
    ppn: 10

Configuration
================================================================================

All configuration of UMD-LETKF is done through a single `YAML <https://yaml.org/spec/1.1/>`_
configuration file. The UMD-LETKF, when run, will by default look for a configuration
file in the same directory named `letkf.yaml`. Or, a different file can be specified
on the command line
::

  ./letkfdriver <somefile.yaml>


An introduction to the YAML format can be found `here <https://yaml.org/spec/1.1/>`_.

.. warning::
   In some places the UMD-LETKF does not perform extensive error checking on the
   configuration file, meaning an improperly defined configuration file might result
   in cryptic messages and crash. If you come across this, please 
   :ref:`let me know <support>` as I am slowly trying to make sure all yaml
   misconfigurations produce helpful error messages.


Each of the major sections of the configuration file are described below. **Each main section is required**, though within a section specific parameters may be optional. 


.. toctree::
   :caption: Configuration Sections:
   :maxdepth: 1

   localization
   mpi
   observation
   solver
   state
   example




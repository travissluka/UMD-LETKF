[![Build Status](https://travis-ci.org/travissluka/UMD-LETKF.svg?branch=develop)](https://travis-ci.org/travissluka/UMD-LETKF)
[![codecov](https://codecov.io/gh/travissluka/UMD-LETKF/branch/develop/graph/badge.svg)](https://codecov.io/gh/travissluka/UMD-LETKF)
[![docs](https://readthedocs.org/projects/umd-letkf/badge/?version=latest)](http://umd-letkf.readthedocs.io)

For complete set of documentation, visit [umd-letkf.readthedocs.io](http://umd-letkf.readthedocs.io/).


```diff
-NOTE: This library is still under active development, though it should have sufficient
-  documentation to be useable now. It has been extensively tested with ocean data assimilation
-  specifically. Appropriate localization classes still need to be added for non-ocean domains.
```

Brief Description
----------
The following library is a rewrite of the local ensemble transform Kalman filter (LETKF) originally developed by Hunt et al., 2007 [1], coded by Takemasa Miyoshi [2], with additional modifications for the ocean by Steve Penny [3].

It is built with the following design choices in mind:

* **model agnostic library** - A single generic LETKF library is provided that can be compiled once and then used in all domains of a coupled LETKF system. Redundancies in code are elimited this way. Most specialization for a given domain are done through configuration files, and a generic driver is provided that should handle most use cases. A custom driver can easily be built to interface with the library if model specific code needs to be added.
* **object oriented design** - Several default implmentations of classes for observation I/O, model state I/O, and localization are provided. If different functionality is required, the user can create their own derived classes and register them with the LETKF library.
* **multi-model strong coupling** - By being model agnostic, the code should allow for easy transition from weakly coupled to strongly coupled DA. The same LETKF code can be used for multiple independent executables (one for each domain), and cross-domain observations can be assimilated by selecting the appropriate observation I/O and localization classes.

[1]: Hunt, B. R., Kostelich, E. J., & Szunyogh, I. (2007). Efficient data assimilation for spatiotemporal chaos: A local ensemble transform Kalman filter. Physica D: Nonlinear Phenomena, 230(1-2), 112–126. [http://doi.org/10.1016/j.physd.2006.11.008](http://doi.org/10.1016/j.physd.2006.11.008)

[2]: Miyoshi, T. (2005). Ensemble Kalman Filter Experiments with a Primitive-Equation Global Model. University of Maryland. Retrieved from [http://hdl.handle.net/1903/3046](http://hdl.handle.net/1903/3046)

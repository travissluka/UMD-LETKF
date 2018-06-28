[![Build Status](https://travis-ci.org/travissluka/UMD-LETKF.svg?branch=develop)](https://travis-ci.org/travissluka/UMD-LETKF)
[![codecov](https://codecov.io/gh/travissluka/UMD-LETKF/branch/develop/graph/badge.svg)](https://codecov.io/gh/travissluka/UMD-LETKF)

Brief Description
----------
The following library is a rewrite of the local ensemble transform Kalman filter (LETKF) originally developed by Hunt et al., 2007 [1], coded by Takemasa Miyoshi [2], with additional modifications by Steve Penny [3]. It should be seen as an intermediate step of generalizing the current LETKF code before incorporating into the JEDI project for NCEP. 

It is built with the following design choices in mind:

* **model agnostic library** -
* **object oriented design** - 
* **multi-model strong coupling** - 

[1]: Hunt, B. R., Kostelich, E. J., & Szunyogh, I. (2007). Efficient data assimilation for spatiotemporal chaos: A local ensemble transform Kalman filter. Physica D: Nonlinear Phenomena, 230(1-2), 112–126. [http://doi.org/10.1016/j.physd.2006.11.008](http://doi.org/10.1016/j.physd.2006.11.008)

[2]: Miyoshi, T. (2005). Ensemble Kalman Filter Experiments with a Primitive-Equation Global Model. University of Maryland. Retrieved from [http://hdl.handle.net/1903/3046](http://hdl.handle.net/1903/3046)

[3]
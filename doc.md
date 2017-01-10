---
project: UMD-LETKF
src_dir: ./src
output_dir: ./doc
project_github: https://github.com/travissluka/umd-letkf
summary: Universal Multi-Domain Local Ensemble Transform Kalman Filter
author: Travis Sluka
email: tsluka@umd.edu
github: https://github.com/travissluka
graph: true
sort: permission
source: false
display: public
display: protected
display: private
extra_mods: netcdf : http://www.unidata.ucar.edu/software/netcdf/docs
page_dir: doc_src/pages
---

Brief Description
----------
The following library is a rewrite of the local ensemble transform Kalman filter (LETKF) originally developed by Hunt et al., 2007 [^1] and coded by Takemasa Miyoshi [^2]. It is built with the following design choices in mind:

* **model independent generalized library** -
* **object oriented design** -
* **multi-model strong coupling** -

[^1]: Hunt, B. R., Kostelich, E. J., & Szunyogh, I. (2007). Efficient data assimilation for spatiotemporal chaos: A local ensemble transform Kalman filter. Physica D: Nonlinear Phenomena, 230(1-2), 112–126. [http://doi.org/10.1016/j.physd.2006.11.008](http://doi.org/10.1016/j.physd.2006.11.008)

[^2]: Miyoshi, T. (2005). Ensemble Kalman Filter Experiments with a Primitive-Equation Global Model. University of Maryland. Retrieved from [http://hdl.handle.net/1903/3046](http://hdl.handle.net/1903/3046)

*[LETKF]: Local Ensemble Transform Kalman Filter
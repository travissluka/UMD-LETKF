Inflation Schemes
=======================

In order to account for an underdispersive ensemble, several multiplicative inflation schemes have been implemented in UMD-LETKF (with hopefully more to be implemented eventually). If you're not sure which one to pick, it is usually safest to choose :ref:`RTPS` with a value between 0.0 and 1.0.



.. _mul_infl:

Multiplicative
-----------------
The inflation factor :math:`\alpha`, which is greater than or equal to 1.0, increases the magnitude of the analysis perturbations.

.. math::
   \mathbf{x}_i^{'a} \leftarrow \alpha \mathbf{x}_i^{'a}

This method works sufficiently for domain that are regulary sampled by observations. (e.g. the atmosphere). If a domain is **not** sufficiently sampled (such as the deep ocean), this method may result in the ensemble spread growing far too rapidly and the filter ultimately diverging.


.. _RTPP:
   
Relaxation to Prior Perturbations (RTPP)
------------------------------------------
The perturbations of the analysis, :math:`\mathbf{x}_i^{'a}` are relaxed a percentage, :math:`\alpha`, back to the background perturbations, :math:`\mathbf{x}_i^{'b}` [Zhang2004]_. This has the benefit of effectively being a combination of both multiplicative, and additive inflation.

.. math::
   \mathbf{x}_i^{'a} \leftarrow \left ( 1 - \alpha \right ) \mathbf{x}_i^{'a} + \alpha \mathbf{x}_i^{'b}

   

.. _RTPS:

Relaxation to Prior Spread (RTPS)
-----------------------------------
The *spread* of the analysis, :math:`\sigma^a`, is relaxed a percentage of the way, :math:`\alpha`, back to the spread of the background, :math:`\sigma^b` [Whitaker2012]_.

.. math::
   \mathbf{x}_i^{'a} \leftarrow \mathbf{x}_i^{'a} \left ( \alpha \frac{\sigma^b - \sigma^a}{\sigma_b} +1 \right )

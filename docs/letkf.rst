LETKF Algorithm
====================

The following is a brief conceptual overview from [Sluka2016]_ of how the LETKF algorithm operates, for a complete description see [Hunt2007]_.

The local ensemble transform Kalman filter (LETKF) is a type of ensemble Kalman filter (EnKF) which uses an ensemble of forecasts :math:`\left\{\mathbf{x}^{b(i)} : i = 1,2,...,k \right\}` to determine the statistics of the background error covariance. This information is combined with new observations, :math:`\mathbf{y}^o`, to generate an analysis mean, :math:`\bar{\mathbf{x}}^a`, and a set of new ensemble members, :math:`\mathbf{x}^{a(i)}`. First, the model state is mapped to observation space by applying a nonlinear observation operator :math:`H` to each background ensemble member

.. math::
   \mathbf{y}^{b(i)} = H\mathbf{x}^{b(i)}

note, that the application of the observation operator is applie *outside* this UMD-LETKF library.

A set of intermediate weights, :math:`\bar{\mathbf{w}}^{a}` are calculated to find the analysis mean :math:`\bar{\mathbf{x}}^a`

.. math::
   \tilde{\mathbf{P}}^a =
   \left [
     \left( k-1 \right ) \mathbf{I} +
     \left( \mathbf{Y}^b \right )^T \mathbf{R}^{-1} \mathbf{Y}^b
   \right ]^{-1}

.. math::
   \bar{\mathbf{w}}^a =
   \tilde{\mathbf{P}}^a \left( \mathbf{Y}^b \right)^T \mathbf{R}^{-1}
   \left( \mathbf{y}^o - \bar{\mathbf{y}}^b \right)

.. math::
   \bar{\mathbf{x}}^a =
   \bar{\mathbf{x}}^b + \mathbf{X}^b \bar{\mathbf{w}}^a

where :math:`\bar{\mathbf{x}}^b` and :math:`\bar{\mathbf{y}}^b` are the ensemble mean of the background in model space and observation space, respectively. :math:`\mathbf{X}^b` and :math:`\mathbf{Y}^b` are the matrices whose columns represent the ensemble perturbations from those means, and :math:`\mathbf{R}` is the observation error covariance matrix.

Last, the set of intermediate weights, :math:`\mathbf{W}^a` are calculated to find the perturbations in model space for the analysis ensemble by

.. math::
   \mathbf{W}^a = \left[ \left( k-1 \right) \tilde{\mathbf{P}}^a \right]^{1/2}

.. math::
   \mathbf{X}^a = \mathbf{X}^b \mathbf{W}^a

the final analysis ensemble members, :math:`\mathbf{x}^{a(i)}`, are the result of adding each column of :math:`\mathbf{X}^a` to :math:`\bar{\mathbf{x}}^a`

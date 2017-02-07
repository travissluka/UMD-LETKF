title: Inflation

Description of inflation types available by default:

How to add own custom inflation type:

### Constant multiplicative ###

\( \mathbf{X}^{'a}_i \leftarrow \alpha \mathbf{X}^{'a}_i \)



### Relaxation to prior spread [^1] (RTPS)  ###
A purely multiplicative inflation that increases the ensemble standard deviation of the analysis some percentage, \( \alpha \), back towards the background ensemble standard deviation.

\( \mathbf{X}^{'a}_i \leftarrow  \mathbf{X}^{'a}_i ( \alpha \cfrac{\sigma^b - \sigma^a}{\sigma^a}  + 1 ) \)


[^1]: Whitaker, J. & T. Hamil, (2015). Evaluating Methods to Account for System Errors in Ensemble Data Assimilation. Monthly Weather Review, 140(9), 3078-3089.[http://doi.org/10.1175/MWR-D-11-00276.1](http://doi.org/10.1175/MWR-D-11-00276.1)


### Relaxation to prior pertubations [^2] (RTPP) ###
Relaxes the analysis perturbations a percentage, \( \alpha \), of the way back toward the background pertubations. This effectively results in a combination of both additive and multiplicative inflation. However, the range of values for \(\alpha\) that work well are likey to be smaller than for RTPS [^1]

\( \mathbf{x}^{'a}_i \leftarrow (1-\alpha) \mathbf{x}^{'a}_i + \alpha \mathbf{x}^{'b}_i \)

[^2]: Zhang, F., C. Snyder, & J. Sun, (2014). Impacts of initial estimate and observation availability on convective-scale data assimilation with an ensemble Kalman filter. Monthly Weather Review, 131, 1238-1253.
# ETKFR

## An Ensemble transform Kalman filter with observation error covariance estimation

This repository includes software for the Ensemble transform Kalman filter with online time-dependent observation error covariance estimation developed in Waller (2013) and Waller et al. (2014). The method allows spatially correlated and time-dependent observation error to be diagnosed and incorporated in an ensemble data assimilation system. The method combines an ensemble transform Kalman filter with a method that uses statistical averages of background and analysis innovations (the so-called Desroziers diagnostic (Desroziers et al., 2005) to provide an estimate of the observation error covariance matrix. 

## References

1. Desroziers, G., Berre, L., Chapnik, B. and Poli, P. (2005) Diagnosis of observation, background and analysis-
error statistics in observation space. Quarterly Journal of the Royal Meteorological Society, 131:
3385–3396, 2005.

1. Kassam, A. and Trefethen, L. (2005) Fourth-order time-stepping for stiff pdes. SIAM J. Sci. Computing,
25:1214–1233.

1. Waller, J. A. (2013) Using observations at different spatial scales in data assimilation for environmental
prediction. PhD thesis, University of Reading, Department of Mathematics and Statistics.
http://www.reading.ac.uk/maths-and-stats/research/theses/maths-phdtheses.aspx.

1. Waller, J. A., Dance, S. L., Lawless, A. S. and Nichols, N. K. (2014) Estimating correlated observation error statistics
 using an ensemble transform Kalman filter. Tellus A, 66. 23294. ISSN 1600-0870 doi: https://doi.org/10.3402/tellusa.v66.23294

## Overview of ETKFR source code 
The ETKFR software is written in Matlab and consits of three functions. The main function is `ETKFR.m` and this in turn calls
the functions `enkf_forecast.m` which forecasts the ensemble members and `etkf_analysis.m` which performs the analysis step. 

We provide an additional mfile `Homogeneous.m` a regularisation method that may be used to regularise the estimated matrix.

We now give specific details of each function.

### ETKFR 

`[dob,doa,MEAN XF,MEAN XA,REst] = ETKFR(model,t,kl,xf0,y,H,AllR,s,regularisation)` . Runs the ETKF assimilation and forecast; it calculates and stores the background and analysis innovations and uses them to calculate the correlated error covariance matrix R.

*Inputs:*
1. `model` - The forecast model. Should have form `x1 = model(t0, t1, x0)`. 
    
    *Inputs:*
    
    1. `t0` - Initial time.
    1. `t1` - Vector of times when model solution required.
    1. `x0` - State vector at initial time.
    
    *Outputs:*
    
    1. `x1` - State vectors at times supplied in t1.
    
1. `t` - A vector of forecast times.
1. `kl` - A vector (length p) of forecast times at which the observations are available.
1. `xf0` - The initial forecast ensemble as a matrix of column vectors (Nm × N).
1. `y` - Array (p × Np) of observation vectors.
1. `H` - Observation operator (Np × Nm).
1. `AllR` - A 3D array (Np ×Np ×Ns) containing the observation error covariance matrix to be used at each step of the assimilation (before the estimated matrix is used).
1. `Ns` - The number of samples that will be used to calculate the observation error covariance
matrix.
1. `regularisation` - The method for regularising the estimated observation error covariance
matrix. Should have the form `[RegR] = regularisation(R)`. 

    *Inputs:*
    
    1. `R` - matrix to be regularised.
    
    *Outputs:*
    
    1. `RegR`- Regularised matrix.
    
*Outputs:*
1. `dob` - The matrix (p × Np) of background innovations.
1. `doa` - The matrix (p × Np) of analysis innovations.
1. `MEAN_XF`- The forecast mean (Nm × p) at each assimilation step.
1. `MEAN_XA`- The analysis mean (Nm × p) at each assimilation step.
1. `EstR` - The estimated error covariance (Nm × p) for assimilation steps after the first Ns
assimilations.

### enkf_forecast
`[xf] = enkf_forecast(model, t0, t1, xa)` Performs the forecast step of the ETKF.

*Inputs:*
1. `model` - The forecast model.
1. `t0` - Initial time.
1. `t1` - Vector of times when model solution required.
1. `xa` - Analysis state vectors to be forecast for each ensemble member.

*Outputs:*
1. `xf` - Forecast state vectors for each ensemble member.

### etfk analysis
`[xa,mean_xf,mean_xa] = etkf_analysis(xf, y, H, rR)`

*Inputs:*
1. `xf` - Forecast state vectors for each ensemble member.
1. `y` - Vector (1 × Np) of observations for the assimilation step
1. `H` - Observation operator, a matrix of size Np × Nm
1. `rR` - Analysis state vectors to be forecast for each ensemble member.

*Outputs:*
1. `xa` - Forecast state vectors to be forecast for each ensemble member.
1. `mean_xf` - Forecast state vectors to be forecast for each ensemble member.
1. `mean_xa` - Forecast state vectors to be forecast for each ensemble member.

### Homogeneous
`[RegR] = Homogeneous(R)`. Code to regularise the estimated observation error matrix by making
it isotropic and homogeneous.
1. *Inputs:*
    1. `R` - matrix to be regularised.
1. *Outputs:*
    1. `RegR` - Regularised matrix.



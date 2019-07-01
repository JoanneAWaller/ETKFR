
# ETKFR

## An Ensemble transform Kalman filter with observation error covariance estimation
### Overview 
This repository includes software for the Ensemble transform Kalman filter with online time-dependent observation error covariance estimation developed in Waller (2013) and Waller et al. (2014). The method allows spatially correlated and time-dependent observation error to be diagnosed and incorporated in an ensemble data assimilation system. The method combines an ensemble transform Kalman filter with a method that uses statistical averages of background and analysis innovations (the so-called Desroziers diagnostic, henceforth denothed the DBCP diagnostic (Desroziers et al., 2005)) to provide an estimate of the observation error covariance matrix. 

### Background

Data assimilation techniques combine physical observations of a dynamical system with a model prediction of the state of the system, known as the background state, weighted by their respective error covariance matrices, to obtain a best estimate of the current state of the system, known as the analysis. We denote here ***x***<sup>*f*</sup> &isin; &reals;<sup>*N*<sup>*m*</sup></sup> as the background state, ***x***<sup>*a*</sup> &isin; &reals;<sup>*N*<sup>*m*</sup></sup> as the analysis and ***y*** &isin; &reals;<sup>*N*<sup>*p*</sup></sup> as the observationsl; the covariance matrices of the errors in the background and observations by ***P***<sup>f</sup> and ***R*** respectively. The background at the next step of the assimilation process is generated by evolving the (possibly nonlinear) dynamical model of the system ***M*** forward from the analysis.

Until recently the observation error covariance matrix has been assumed uncorrelated. However, it has been shown that the error is correlated Waller et al. (2013) and that the inclusion of the correlated errors in the assimilation leads to: a more accurate analysis, the inclusion of more observation information content and an improvement in the NWP skill score Stewart et al. (2013),
Weston et al. (2014).

### Desroziers et al. Diagnostic

The DBCP diagnostic described in Desroziers et al. (2005) shows that the observation error covari-
ance matrix can be calculated using, 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ***R*** &asymp; E[(***y***-***Hx***<sup>a</sup>)(***y***-***Hx***<sup>f</sup>)<sup>T</sup>], (1) 

where the observation operator, ***H***, maps from model to observation space. Equation (1) is valid if the observation and forecast errors used in the gain matrix,

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ***K*** = ***P***<sup>f</sup>***H***<sup>T</sup>(***H******P***<sup>f</sup>***H***<sup>T</sup> + ***R***)<sup>-1</sup>, (2) 

to calculate the analysis, are correct. However, the estimated observation error statistics  may be imporved with successive iterations of the diagnostic (Desroziers et al., 2005). Furthermore, the diagnostic can still provide useful information about the true observation uncertainties (Waller et al., 2016).

### Ensemble transform Kalman filter with R estimation

The ETKFR uses the ETKF and the DBCP diagnostic to estimate a possibly non-uniform time varying observation error covariance matrix. After the filter is initialised the filter is split into two stages. The first is a spin-up stage that runs for a predetermined number of steps, *N*<sup>s</sup>, and is an application of the standard ensemble transform Kalman filter. In the second stage at each assimilation step the observation error covariance matrix is updated using the DBCP diagnostic. We now present in detail the method that we have developed. Here the observation operator, **H**, is chosen to be linear, but the method could be extended to account for a non-linear observation operator ***H*** (e.g. Evensen (2003)). 

The steps of the ETKFR are as follows (note, equations are given in Waller et al. (2014)):


**Initialisation** - Begin with an initial ensemble ***x***<sub>0</sub><sup>i</sup> for *i* = 1...*N* at time *t* = 0 that has an associated initial covariance matrix ***P***<sup>*f*</sup>. Also assume an initial estimate of the observation error covariance matrix ***R***<sub>0</sub>; it is possible that this could just consist of the instrument error.

**Step 1** - Use the full non-linear model to forecast each ensemble member.

**Step 2** - Calculate the forecast ensemble mean and covariance.

**Step 3** - Using the ensemble mean and the observations at time *t*<sub>n</sub>, calculate and store the observation-minus-background departures, ***d***<sup>b</sup> = ***y*** − **H*****x***<sup>*f*</sup>. 

**Step 4** - Update the ensemble mean.

**Step 5** - Calculate the analysis perturbations. 

**Step 6** - Use the analysis mean to calculate the sobservation-minus-analysis departures, ***d***<sup>a</sup> = ***y*** − **H*****x***<sup>*a*</sup>. 

**Step 7** - If *n* > *N*<sup>s</sup>, where *N*<sup>s</sup> is the specified sample size, update ***R*** using the observation-minus-analyis and observation-minus-background departures from the previous *N*<sup>s</sup> assimialtion steps, then symmetrise the matrix. Otherwise keep ***R***<sub>n+1</sub> = ***R***<sub>0</sub>.

Many of the steps in the proposed method are identical to the ETKF. Step 7, along with the storage of the background and analysis innovations in steps 3 and 6, are the additions to the ETKF that provide the estimate of the observation error covariance matrix.

In practice the number of samples available will be limited and therefore the estimated observation error covariance matrix will not be full rank. In this case it may be necessary to apply some form of regularisation to the estimated matrix.

### References

1. Desroziers, G., Berre, L., Chapnik, B. and Poli, P. (2005) Diagnosis of observation, background and analysis-
error statistics in observation space. Quarterly Journal of the Royal Meteorological Society, 131:
3385–3396, 2005.

1. Evensen, G. (2003) The ensemble Kalman filter: Theoretical formulation and practical implementation.
Ocean Dynamics, 53:343–367.

1. Stewart, L. M., Dance, S. L. and Nichols, N. K. (2013) Data assimilation with correlated observation
errors: experiments with a 1-D shallow water model. Tellus A, 65.

1. Waller, J. A. (2013) Using observations at different spatial scales in data assimilation for environmental
prediction. PhD thesis, University of Reading, Department of Mathematics and Statistics.
http://www.reading.ac.uk/maths-and-stats/research/theses/maths-phdtheses.aspx.

1. Waller, J. A., Dance, S. L., Lawless, A. S. and Nichols, N. K. (2014) Estimating correlated observation error statistics
 using an ensemble transform Kalman filter. Tellus A, 66. 23294. ISSN 1600-0870 doi: https://doi.org/10.3402/tellusa.v66.23294
 
1. Waller, J. A., Dance, S. L. and Nichols, N. K. (2016) Theoretical insight into diagnosing observation error correlations using observation-minus-background and observation-minus-analysis statistics. Quarterly Journal of the Royal Meteorological Society, 142 (694). pp. 418-431. ISSN 1477-870X doi: https://doi.org/10.1002/qj.2661
 
1. Weston, P. P. , Bell, W. and Eyre, J. R. (2014) Accounting for correlated error in the assimilaiton of high
resolution sounder data. Quarterly Journal of the Royal Meteorological Society. DOI: 10.1002/qj.2306.

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

*Inputs:*
1. `R` - matrix to be regularised.

*Outputs:*
1. `RegR` - Regularised matrix.



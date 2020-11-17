It contains R code for WN test of high-dimensional time series with observations X: nxp (n realizations of p-dimensional time series). It can also do the WN test for fitted residuals of VAR model based on LASSO estimation.  The main code is wnSVAR.R

Example

    source('wnSVAR.R')
    X = matrix(rnorm(100*20), 100, 20)
    output = wnSVAR(X, K = 10)
    
    output$pvalues1  # p-values at lags 1--K, based on infinity-norm (maxi-norm)
    output$pvalues2  # p-values at lags 1--K, based on L2-norm (maxi-norm)
    
    output$pvalues1   # p-values at lags 1--K, based on infinity-norm (maxi-norm) for sn = p
    output$pvalues10  # p-values at lags 1--K, based on L2-norm (maxi-norm) for sn = p
    
    
    
    



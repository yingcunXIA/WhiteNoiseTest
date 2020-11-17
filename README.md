It contains R code for WN test of high-dimensional time series with observations X: nxp (n realizations of p-dimensional time series). It can also do the WN test for fitted residuals of VAR model based on LASSO estimation.  The main code is wnSVAR.R

It also contans a code for Chernozhukov, Chetverikov, and Kato (2019), ``Inference on causal and structural parameters using many moment inequalities,”
Review of Economic Studies. The main code is CCK.R


Example

    source('wnSVAR.R')  # our proposed method
    X = matrix(rnorm(100*20), 100, 20)
    output = wnSVAR(X, K = 10)
    
    output$pValue0  # p-values at lags 1--K, based on infinity-norm (maxi-norm)
    output$pValue2  # p-values at lags 1--K, based on L2-norm (maxi-norm)
    
    output$pValue0sn   # p-values at lags 1--K, based on infinity-norm (maxi-norm) for sn = p
    output$pValue2sn  # p-values at lags 1--K, based on L2-norm (maxi-norm) for sn = p
    
    source('CCK.R')  # based on Chernozhukov et al (2019)
    out.cck = CCK(X)
    out.cck$pValue0  # p-values at lags 1--K, based on infinity-norm (maxi-norm)
    out.cck$pValue2  # p-values at lags 1--K, based on L2-norm (maxi-norm)
    
    
    



CCK <- function(X, K= 10, Boot=1250, rsample="B")
{

#  source('Macf.R')

  p = dim(X)[2]
  for (i in 1:p)
    X[,i] = (X[,i]-mean(X[,i]))/(sd(X[,i])+1.0e-10)

    np = dim(X);
    n = np[1]
    p = np[2]

    T2 = rep(0, K)
    T0 = rep(0, K) 
    Ts1 = rep(0, K)
    Ts10 = rep(0, K)
    
    for (i in 1:K)
    {
      I1 = 1:(n-i)
      I2 = (1+i):n
      s = t(X[I1,])%*%X[I2,] /sqrt(n)
      T2[i] =  sum(abs(s))
      T0[i] =  max(abs(s))
    }
    
    T2 = cumsum(T2)
    T0 = cumsum(T0)

#    best combination
#    q = floor(sqrt(n)/2)
#    r = floor(sqrt(n)/2)
#    -1, 1   sampling  
        
    q = floor(sqrt(n)/2)
    r = floor(sqrt(n)/2)
    m = floor(n/(q+r))
    

    S = array(0, dim=c(p, p, m, K))
    for (k in 1:K)
    {
      mn = 0
      for (i in 1:m)
      {
        I = ((i-1)*(q+r)+1):((i-1)*(q+r)+q);
        if (max(I)+m <= n)
        {
          S[,, i, k] = t(X[I,])%*%X[I+k,]/sqrt(q);
          mn = mn +1
        }
      }
#      S[,, i, k] = S[,, i, k]/sqrt(max(mn,1))
    }
    
    Tb2 = matrix(0, Boot, K)
    Tb0 = matrix(0, Boot, K)
    for (ib in 1:Boot)
    {
      if (rsample=="B")
      {
        e = (rnorm(m)>0)*2-1
      }
      else
      {
        e = rnorm(m)
      }
      
      for (k in 1:K)
      {
        M = 0
        for (i in 1:m)
        {
          M = M + e[i]*S[,,i, k]/sqrt(m); 
        }
        Tb2[ib,k] = sum(abs(M))
        Tb0[ib,k] = max(abs(M))
      }
      Tb2[ib,] = cumsum(Tb2[ib,])
      Tb0[ib,] = cumsum(Tb0[ib,])
    }
    

#  pValue2 = mean(Tb2>T2);
#  pValue0 = mean(Tb0>T0);
    
    pValue2 = rep(0, K)
    pValue0 = rep(0, K)
    for (k in 1:K)
    {
      pValue2[k] = mean(Tb2[,k]>T2[k])
      pValue0[k] = mean(Tb0[,k]>T0[k])
    }

  return(list(pValue2=pValue2, pValue0=pValue0))
}

Macf <- function(y,K)
{
  np = dim(y);
  n = np[1]
  p = np[2]
  
  
  y = y - t(matrix(rep(colMeans(y), n), p, n));
  
  ACF2 = c();
  ACF0 = c();
  ACFs1 = c();
  ACFs10 = c();
  
  S = 1/sqrt(colSums(y^2))
  if (p > 1) S = diag(S, p)
  
  if ((p > 1)& (p<n/5))    ######################
  {
    S = t(y)%*%y;
    S = svd(S)
    S = S$u%*%diag(1/(sqrt(S$d)+1/n^3))%*%t(S$u)
  }
  
  acf0 = 0;
  acf2 = 0;
  
  acf_s1 = 0;
  acf_s10 = 0;
  
  rA = 0
  #  sn = min(floor(n/log(p)^2), p^2)
  sn = p
  
  for (i in 1:K)
  {
    rk = abs(S%*% (t(y[1:(n-i), ])%*%y[(i+1):n, ]) %*% S)/(n-i)
    rA = rA + rk
    si = rowSums(rk)
    
    srk = sort(rA, decreasing = TRUE) 
    acf0 = sum(srk[1]);    
    acf2 = acf2 + sum(si);    
    
    
    rk = sort(rk, decreasing = TRUE)
    acf_s1 =  max(acf_s1, sum(srk[1:sn]));
    acf_s10 = acf_s10 + sum(rk[1:sn]);
    
    
    ACF2 = c(ACF2, acf2);
    ACF0 = c(ACF0, acf0);
    ACFs1 = c(ACFs1, acf_s1);
    ACFs10 = c(ACFs10, acf_s10);
  }
  return(list(ACF2=ACF2, ACF0=ACF0, ACFs1=ACFs1, ACFs10=ACFs10))
}


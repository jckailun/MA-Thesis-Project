##### Newton's method for cohorts #####
## 0. Init. Parameters 
G = 3; cohort = 1
dat = clust(cohort, G)
df = dat[[1]]
G = length(df) # = G
D = dat[[2]]
theta0 = rep(0, 3*k + G)

## 1. Log-lik
ll = function(theta, G, k = 16, cohort, BS = FALSE){ # validated
  
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]
  l = 0
  
  for (j in 1:G){
    
    d = as.matrix(df[[j]])
    aj = a[j]
    n = nrow(d)
    D1 = D[[j]][,1]; D2 = D[[j]][,2]; D3 = D[[j]][,3]; D4 = D[[j]][,4]   
    
    for (i in 1:n){
      
      x_i = as.vector(d[i,])
      eta1 = as.numeric(x_i %*% b1 + aj)
      eta2 = as.numeric(x_i %*% b2 + aj)
      eta3 = as.numeric(x_i %*% b3 + aj)
      
      t1 = -log(1 + exp(eta1))*D1[i]
      t2 = (eta1 - log(1 + exp(eta1)) - log(1 + exp(eta2))) * D2[i]
      t3 = (eta1 + eta2 - log(1 + exp(eta1)) - log(1 + exp(eta2)) - log(1 + exp(eta3))) * D3[i]
      t4 = (eta1 + eta2 + eta3 - log(1 + exp(eta1)) - log(1 + exp(eta2)) - log(1 + exp(eta3))) * D4[i]
      
      l = l + (t1 + t2 + t3 + t4)
      
    }
  }
  l
}
ll(theta0, G, k=16, cohort=1)

## 2. Score Function
dldb1 = function(theta, G, k = 16, cohort){ # validated

  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]

  B1 = 0

  for (j in 1:G){
    
    d = as.matrix(df[[j]])
    aj = a[j]
    n = nrow(d)
    
    D1 = D[[j]][,1]; D2 = D[[j]][,2]; D3 = D[[j]][,3]; D4 = D[[j]][,4]   
    
    for (i in 1:n){
      x_i = as.vector(d[i,])
      eta1 = as.numeric(x_i %*% b1 + aj)
      #eta2 = as.numeric(x_i %*% b2 + aj)
      #eta3 = as.numeric(x_i %*% b3 + aj)
      
      t1 = -plogis(eta1) * x_i * D1[i]
      t2 = (1-plogis(eta1)) * x_i * D2[i]
      t3 = (1-plogis(eta1)) * x_i * D3[i]
      t4 = (1-plogis(eta1)) * x_i * D4[i]
      
      B1 = B1 + (t1 + t2 + t3 + t4)
    }
  }
  
  as.vector(B1)
  
}
dldb1(theta0, G, k=16, cohort=1)

dldb2 = function(theta, G, k=16, cohort){ # validated
  
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]
  
  B2 = 0
  
  for (j in 1:G){
    
    d = as.matrix(df[[j]])
    aj = a[j]
    n = nrow(d)
    
    D1 = D[[j]][,1]; D2 = D[[j]][,2]; D3 = D[[j]][,3]; D4 = D[[j]][,4]
    
    for (i in 1:n){
      x_i = as.vector(d[i,])
      #eta1 = as.numeric(x_i %*% b1 + aj)
      eta2 = as.numeric(x_i %*% b2 + aj)
      #eta3 = as.numeric(x_i %*% b3 + aj)
      
      t2 = -plogis(eta2) * x_i * D2[i]
      t3 = (1-plogis(eta2)) * x_i * D3[i]
      t4 = (1-plogis(eta2)) * x_i * D4[i]
      
      B2 = B2 + (t2 + t3 + t4)
    }
  }
  
  as.vector(B2)
  
}
dldb2(theta0, G, k = 16, cohort = 1)

dldb3 = function(theta, G, k=16, cohort){ # validated
  
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]
  
  B3 = 0
  
  for (j in 1:G){
    
    d = as.matrix(df[[j]])
    aj = a[j]
    n = nrow(d)
    D1 = D[[j]][,1]; D2 = D[[j]][,2]; D3 = D[[j]][,3]; D4 = D[[j]][,4]
    
    for (i in 1:n){
      x_i = as.vector(d[i,])
      #eta1 = as.numeric(x_i %*% b1 + aj)
      #eta2 = as.numeric(x_i %*% b2 + aj)
      eta3 = as.numeric(x_i %*% b3 + aj)
      
      t3 = -plogis(eta3) * x_i * D3[i]
      t4 = (1-plogis(eta3)) * x_i * D4[i]
      
      B3 = B3 + (t3 + t4)
    }
  }
  
  as.vector(B3)
  
}
dldb3(theta0, G, k=16, cohort=1)

dldag = function(theta, j, k = 16, cohort){ # validated
  
  # g controls groups
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]
  aj = a[j]
  Ag = 0

  D1 = D[[j]][,1]; D2 = D[[j]][,2]; D3 = D[[j]][,3]; D4 = D[[j]][,4]
  
  d = as.matrix(df[[j]])
  n = nrow(d)
  
  for (i in 1:n){
    x_i = as.vector(d[i,])
    eta1 = as.numeric(x_i %*% b1 + aj)
    eta2 = as.numeric(x_i %*% b2 + aj)
    eta3 = as.numeric(x_i %*% b3 + aj)
    
    t1 = - plogis(eta1) * D1[i]
    t2 = (1 - plogis(eta1) - plogis(eta2))  * D2[i]
    t3 = (2 - plogis(eta1) - plogis(eta2) - plogis(eta3)) * D3[i]
    t4 = (3 - plogis(eta1) - plogis(eta2) - plogis(eta3)) * D4[i]
    
    Ag = Ag + (t1 + t2 + t3 + t4)
  }
  Ag
}

dldag(theta0, 1, k=16, cohort=1)
dldag(theta0, 2, k=16, cohort=1)
dldag(theta0, 3, k=16, cohort=1)

score = function(theta, G, k=16, cohort){ # validated
  
  sa = NULL
  
  sb1 = dldb1(theta, G, k = 16, cohort)
  sb2 = dldb2(theta, G, k = 16, cohort)
  sb3 = dldb3(theta, G, k = 16, cohort)
  
  for (j in 1:G){
    sa = c(sa, dldag(theta, j, k = 16, cohort))
  }
  
  c(sb1, sb2, sb3, sa)
  
}
score(theta0, G, k=16, cohort=1)

## 3. Hessian Matrix ##
d2ldb12 = function(theta, G, k = 16, cohort){ # validated
  
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]
  BB1 = 0

  for (j in 1:G){
    
    d = as.matrix(df[[j]])
    ag = a[j]
    n = nrow(d)
    
    for (i in 1:n){
      
      x_i = as.vector(d[i,])
      eta1 = as.numeric(x_i %*% b1 + ag)
      t1 = plogis(eta1) * (1-plogis(eta1)) * (x_i%o%x_i)
      BB1 = BB1 + t1
      
    }
    
  }
  
  -BB1
  
}
d2ldb12(theta0, G, k=16, 1)

d2ldb22 = function(theta, G, k = 16, cohort){ # validated
  
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]
  BB2 = 0

  for (j in 1:G){
    
    d = as.matrix(df[[j]])
    ag = a[j]
    n = nrow(d)
    D1 = D[[j]][,1]; D2 = D[[j]][,2]; D3 = D[[j]][,3]; D4 = D[[j]][,4]
    
    for (i in 1:n){
      
      x_i = as.vector(d[i,])
      #eta1 = as.numeric(x_i %*% b1 + aj)
      eta2 = as.numeric(x_i %*% b2 + ag)
      #eta3 = as.numeric(x_i %*% b3 + aj)
      
      t1 = (1 - D1[i]) * plogis(eta2) * (1-plogis(eta2)) * (x_i%o%x_i)
      BB2 = BB2 + t1
      
    }
    
  }
  
  -BB2
  
}
d2ldb22(theta0, G, k=16, cohort = 1)

d2ldb32 = function(theta, G, k = 16, cohort){ # validated
  
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  a = theta[-(1:(3*k))]
  BB3 = 0
  
  for (j in 1:G){
    
    d = as.matrix(df[[j]])
    ag = a[j]
    n = nrow(d)
    
    D1 = D[[j]][,1]; D2 = D[[j]][,2]; D3 = D[[j]][,3]; D4 = D[[j]][,4]
    
    for (i in 1:n){
      
      x_i = as.vector(d[i,])
      #eta1 = as.numeric(x_i %*% b1 + aj)
      #eta2 = as.numeric(x_i %*% b2 + aj)
      eta3 = as.numeric(x_i %*% b3 + ag)
      
      t1 = (D3[i] + D4[i]) * plogis(eta3) * (1-plogis(eta3)) * (x_i%o%x_i)
      BB3 = BB3 + t1
      
    }
    
  }
  
  -BB3
  
}
d2ldb32(theta0, G, k=16, cohort)

d2ldag2 = function(theta, j, k = 16, cohort){ # validated
  # g controls groups
  
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  ag = theta[3*k + j]
  
  AAg = 0
  
  d = as.matrix(df[[j]])
  n = nrow(d)
  
  D1 = D[[j]][,1]; D2 = D[[j]][,2]; D3 = D[[j]][,3]; D4 = D[[j]][,4]
  
  for (i in 1:n){
    x_i = as.vector(d[i,])
    eta1 = as.numeric(x_i %*% b1 + ag)
    eta2 = as.numeric(x_i %*% b2 + ag)
    eta3 = as.numeric(x_i %*% b3 + ag)
    
    t1 = plogis(eta1) * (1-plogis(eta1)) * D1[i]
    t2 = (plogis(eta1)*(1-plogis(eta1)) + plogis(eta2)*(1-plogis(eta2))) * D2[i]
    t3 = ((plogis(eta1)*(1-plogis(eta1)) + plogis(eta2)*(1-plogis(eta2)) + plogis(eta3)*(1-plogis(eta3)))) * (D3[i] + D4[i])
    
    AAg = AAg + (t1 + t2 + t3)
  }
  -AAg
}
d2ldag2(theta0, 1, k=16, cohort=1)
d2ldag2(theta0, 2, k=16, cohort=1)
d2ldag2(theta0, 3, k=16, cohort=1)

d2ldb1dag = function(theta, j, k = 16, cohort){ # validated

  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  ag = theta[3*k + j]

  B1A = rep(0, k)
  
  d = as.matrix(df[[j]])
  n = nrow(d)
  
  D1 = D[[j]][,1]; D2 = D[[j]][,2]; D3 = D[[j]][,3]; D4 = D[[j]][,4]

  for (i in 1:n){
    x_i = as.vector(d[i,])
    eta1 = as.numeric(x_i %*% b1 + ag)
    #eta2 = as.numeric(x_i %*% b2 + ag)
    #eta3 = as.numeric(x_i %*% b3 + ag)
    
    t1 = plogis(eta1) * (1-plogis(eta1)) * x_i
    
    B1A = B1A + t1
  }
  -as.matrix(B1A, nr=k, nc=1)
}
d2ldb1dag(theta0, 1, k = 16, cohort)
d2ldb1dag(theta0, 2, k = 16, cohort)
d2ldb1dag(theta0, 3, k = 16, cohort)

d2ldb2dag = function(theta, j, k = 16, cohort){ # validated
  
  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  ag = theta[3*k + j]
  
  G = length(df)
  B2A = rep(0, k)
  d = as.matrix(df[[j]])
  n = nrow(d)
  
  D1 = D[[j]][,1]; D2 = D[[j]][,2]; D3 = D[[j]][,3]; D4 = D[[j]][,4]
  
  for (i in 1:n){
    x_i = as.vector(d[i,])
    #eta1 = as.numeric(x_i %*% b1 + ag)
    eta2 = as.numeric(x_i %*% b2 + ag)
    #eta3 = as.numeric(x_i %*% b3 + ag)
    
    t1 = plogis(eta2) * (1-plogis(eta2)) * x_i * (D2[i] + D3[i] + D4[i])
    
    B2A = B2A + t1
  }
  -as.matrix(B2A, nr=k, nc=1)
}
d2ldb2dag(theta0, 1, k = 16, cohort)
d2ldb2dag(theta0, 2, k = 16, cohort)
d2ldb2dag(theta0, 3, k = 16, cohort)

d2ldb3dag = function(theta, j, k = 16, cohort){ # validated

  b1 = theta[1:k]
  b2 = theta[(k+1):(2*k)]
  b3 = theta[(2*k+1):(3*k)]
  ag = theta[3*k + j]
  
  B3A = rep(0, k)
  d = as.matrix(df[[j]])
  n = nrow(d)
  
  D1 = D[[j]][,1]; D2 = D[[j]][,2]; D3 = D[[j]][,3]; D4 = D[[j]][,4]
  
  for (i in 1:n){
    x_i = as.vector(d[i,])
    #eta1 = as.numeric(x_i %*% b1 + ag)
    #eta2 = as.numeric(x_i %*% b2 + ag)
    eta3 = as.numeric(x_i %*% b3 + ag)
    
    t1 = plogis(eta3) * (1-plogis(eta3)) * x_i * (D3[i] + D4[i])
    
    B3A = B3A + t1
  }
  -as.matrix(B3A, nr=k, nc=1)
}
d2ldb3dag(theta0, j=2, k = 16, cohort=1)

d2ldagdb1 = function(theta, j, k = 16, cohort){ t(d2ldb1dag(theta, j, k = 16, cohort)) } # validated
d2ldagdb2 = function(theta, j, k = 16, cohort){ t(d2ldb2dag(theta, j, k = 16, cohort)) } # validated
d2ldagdb3 = function(theta, j, k = 16, cohort){ t(d2ldb3dag(theta, j, k = 16, cohort)) } # validated

hess = function(theta, G, k = 16, cohort){ # validated
  
  H = matrix(0, length(theta), length(theta))
  H[1:k,1:k] = d2ldb12(theta, G, k, cohort)
  H[((k+1):(2*k)),((k+1):(2*k))] = d2ldb22(theta, G, k, cohort)
  H[((2*k+1):(3*k)),((2*k+1):(3*k))] = d2ldb32(theta, G, k, cohort)
  
  for (j in 1:G){
    H[(3*k)+j, (3*k)+j] = d2ldag2(theta, j, k, cohort)
    
    H[1:k, (3*k)+j] = d2ldb1dag(theta, j, k, cohort)
    H[(k+1):(2*k), (3*k)+j] = d2ldb2dag(theta, j, k, cohort)
    H[(2*k+1):(3*k), (3*k)+j] = d2ldb3dag(theta, j, k, cohort)
    
    H[(3*k)+j, 1:k] = d2ldagdb1(theta, j, k, cohort)
    H[(3*k)+j, (k+1):(2*k)] = d2ldagdb2(theta, j, k, cohort)
    H[(3*k)+j, (2*k+1):(3*k)] = d2ldagdb3(theta, j, k, cohort)
  }
  H
  
}
hess(theta0, G, k=16, cohort)
### End: Hessian Matrix ###

##### Implementation #####
require(pracma)
brd = function(theta, G, k = 16, cohort, maxiter = 100, tol = .Machine$double.eps^(1/2)) {

  x0 = theta
  loglik = ll(x0, G, k, cohort)
  print(loglik)
  
  y0 = score(x0, G, k, cohort)
  A0 = hess(x0, G, k, cohort)
  B0 = ginv(A0)
  
  xnew <- x0 - B0 %*% y0
  ynew <- score(xnew, G, k, cohort)
  niter <- 1
  
  while (niter < maxiter) {
    
    loglik = c(loglik, ll(xnew, G, k, cohort))
    print(ll(xnew, G, k, cohort))
    
    s <- xnew - x0
    d <- ynew - y0
    if (norm(s, "F") < tol || norm(as.matrix(ynew), "F") < tol || is.nan(ll(xnew, G, k, cohort)) == TRUE) 
      break 
    
    print(c(norm(s, "F"), norm(as.matrix(ynew), "F")))
    B0 <- B0 + (s - B0 %*% d) %*% t(s) %*% B0/c(t(s) %*% B0 %*% d)
    x0 <- xnew
    xnew <- xnew - B0 %*% ynew
    y0 <- ynew
    ynew <- score(xnew, G, k, cohort)
    niter <- niter + 1
    
  }
  
  if (niter >= maxiter) 
    warning(paste("Not converged: Max number of iterations reached."))
  fnew <- sqrt(sum(ynew^2))
  return(list(zero = c(xnew), fnorm = fnew, niter = niter, loglik = loglik, hess = hess(xnew, G, k, cohort)))

}




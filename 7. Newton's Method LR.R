##### Newton's method for cohorts #####
theta.lr= rep(0, 48+3)
theta.lr = rep(0, 48+7)

dat = clust(cohort, G, BS = TRUE)

G = 3; cohort = 1
dat = clust(cohort, G)
df = dat[[1]]
G = length(df) # = G
D = dat[[2]]
theta0 = rep(0, 3*k + G)
a = rep(0, G)


ll.lr = function(theta, cohort){ # validated
  
  l = 0
  G = length(theta)

  for (j in 1:G){
    
    aj = theta[j]
    d = as.matrix(df[[j]])
    n = nrow(d)
    D1 = D[[j]][,1]; D2 = D[[j]][,2]; D3 = D[[j]][,3]; D4 = D[[j]][,4]
    
    for (i in 1:n){
      
      eta1 = as.numeric(aj)
      eta2 = as.numeric(aj)
      eta3 = as.numeric(aj)
      
      t1 = -log(1 + exp(eta1)) * D1[i]
      t2 = (eta1 - log(1 + exp(eta1)) - log(1 + exp(eta2))) * D2[i]
      t3 = (eta1 + eta2 - log(1 + exp(eta1)) - log(1 + exp(eta2)) - log(1 + exp(eta3))) * D3[i]
      t4 = (eta1 + eta2 + eta3 - log(1 + exp(eta1)) - log(1 + exp(eta2)) - log(1 + exp(eta3))) * D4[i]
      
      l = l + (t1 + t2 + t3 + t4)
      
    }
    
  }
  
  l
  
}
ll.lr(a, cohort=1)
### Score Function ###

dldag.lr = function(theta, g, cohort){ # validated
  # g controls groups
  
  l = 0
  G = length(theta)
  
  d = as.matrix(df[[g]])
  n = nrow(d)
  D1 = D[[g]][,1]; D2 = D[[g]][,2]; D3 = D[[g]][,3]; D4 = D[[g]][,4]
  
  Ag = 0
  aj = theta[g]
  
  for (i in 1:n){
    
    eta1 = as.numeric(aj)
    eta2 = as.numeric(aj)
    eta3 = as.numeric(aj)
    
    t1 = - plogis(eta1) * D1[i]
    t2 = (1 - plogis(eta1) - plogis(eta2))  * D2[i]
    t3 = (2 - plogis(eta1) - plogis(eta2) - plogis(eta3)) * D3[i]
    t4 = (3 - plogis(eta1) - plogis(eta2) - plogis(eta3)) * D4[i]
    
    Ag = Ag + (t1 + t2 + t3 + t4)
  }
  Ag
}

dldag.lr(a, g=1, cohort=1)

score.lr = function(theta, cohort){ # validated
  
  G = length(theta)
  sa = NULL
  for (j in 1:G){
    sa = c(sa, dldag.lr(theta, g = j, cohort))
  }
  
  sa
  
}
score.lr(a, cohort = 1)
### End: Score Function ###
### Hessian Matrix ###

d2ldag2.lr = function(theta, g, cohort){ # validated

  l = 0
  G = length(theta)
  
  d = as.matrix(df[[g]])
  n = nrow(d)
  D1 = D[[g]][,1]; D2 = D[[g]][,2]; D3 = D[[g]][,3]; D4 = D[[g]][,4]
  
  aj = theta[g]
  
  AAg = 0
  
  for (i in 1:n){
    
    eta1 = as.numeric(aj)
    eta2 = as.numeric(aj)
    eta3 = as.numeric(aj)
    
    t1 = plogis(eta1) * (1-plogis(eta1)) * D1[i]
    t2 = (plogis(eta1) * (1-plogis(eta1)) + plogis(eta2)*(1-plogis(eta2))) * D2[i]
    t3 = ((plogis(eta1) * (1-plogis(eta1)) + plogis(eta2)*(1-plogis(eta2)) + plogis(eta3)*(1-plogis(eta3)))) * (D3[i] + D4[i])
    
    AAg = AAg + (t1 + t2 + t3)
  }
  -AAg
}
d2ldag2.lr(a, g=1, cohort=1)

hess.lr = function(theta, cohort){ # validated
  
  G = length(theta)
  H = matrix(0, length(theta), length(theta))
  
  for (g in 1:G){
    H[g,g] = d2ldag2.lr(theta, g, cohort)
  }
  H
}
hess.lr(a, cohort = 1)
### End: Hessian Matrix ###

##### Implementation #####
require(pracma)
brd.lr = function(theta, cohort, maxiter = 100, tol = .Machine$double.eps^(1/2)) {
  # input df shall be a well-prepared data set
  # standardize whatever it's needed for ESTIMATION
  x0 = theta
  loglik = ll.lr(x0, cohort)
  print(loglik)
  
  y0 = score.lr(x0, cohort)
  A0 = hess.lr(x0, cohort)
  B0 = inv(A0)
  
  xnew <- x0 - B0 %*% y0
  ynew <- score.lr(xnew, cohort)
  niter <- 1
  while (niter < maxiter) {
    
    loglik = c(loglik, ll.lr(xnew, cohort))
    print(ll.lr(xnew, cohort))
    
    s <- xnew - x0
    d <- ynew - y0
    if (norm(s, "F") < tol || norm(as.matrix(ynew), "F") < tol) 
      break
    
    print(c(norm(s, "F"), norm(as.matrix(ynew), "F")))
    B0 <- B0 + (s - B0 %*% d) %*% t(s) %*% B0/c(t(s) %*% B0 %*% d)
    x0 <- xnew
    xnew <- xnew - B0 %*% ynew
    y0 <- ynew
    ynew <- score.lr(xnew, cohort)
    niter <- niter + 1
  }
  if (niter >= maxiter) 
    warning(paste("Not converged: Max number of iterations reached."))
  fnew <- sqrt(sum(ynew^2))
  return(list(zero = c(xnew), fnorm = fnew, niter = niter, loglik = loglik, hess = hess.lr(xnew, cohort)))
}
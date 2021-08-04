##### Estimation with Hyperparameter G > 1 #####
require(MASS)

set.seed(1000)


k = 16
H = c(2, 7)
n = c(nrow(X1.c), nrow(X2.c))
B = 1000

accuracy = as.data.frame(cbind(H, rep(0, length(H)), rep(0, length(H))))
names(accuracy) = c("G", "Model Precision", "Equivalent Bootstrap Size")
out.table = list(matrix(0, nr = 3*k+H[1], nc = 5), matrix(0, nr = 3*k+H[2], nc = 5))
confusion = rep(list(matrix(0, nr = 4, nc = 4)), 2)

for (cohort in 1:2){
  
  X = list(X1.c, X2.c)
  X.c = list(X1.cc, X2.cc)
  diss = list(diss1, diss2)
  dummy = list(dummy.1, dummy.2)
  
  G = H[cohort]
  dat = clust(cohort, G, BS = FALSE)
  df = dat[[1]]
  D = dat[[2]]
  theta0 = rep(0, 3*k + G)
  membership = dat[[3]]

  est01 = brd(theta0, G, k = 16, cohort = cohort)
  est01.bs = NULL
  
  for (b in 1:B){
    
    idx.bs = sample(1:n[cohort], n[cohort], replace = TRUE)
    
    dummy.bs = dummy
    dummy.bs[[cohort]] = dummy.bs[[cohort]][idx.bs, ]
    
    theta.bs = rep(0, 3*k + G)
    
    X.bs = X
    X.bs[[cohort]] = X.bs[[cohort]][idx.bs, ]
    
    X.bs.c = X.c
    X.bs.c[[cohort]] = X.bs.c[[cohort]][idx.bs, ]
    
    gower_dist <- daisy(X.bs.c[[cohort]], metric = "gower", type = list(symm = c(12, 14:17)))
    diss.bs = as.matrix(gower_dist)
    
    dat = clust(cohort, G, BS = TRUE)

    df = dat[[1]]
    D = dat[[2]]

    out = brd(theta.bs, G, k = 16, cohort = cohort)
    est01.bs = rbind(est01.bs, out$zero)
    
  }
  
  mid = apply(est01.bs, 2, median)
  mn = apply(est01.bs, 2, mean)
  
  for (i in 1:B){
    if (abs(mean(est01.bs[i, ])) > max(mn)) { est01.bs[i, ] = NA }
  }
  
  se.bs = sqrt(diag(cov(na.omit(est01.bs))))
  z01 = est01$zero / se.bs

  p01.bs = NULL
  for (v in 1:length(theta0)){
    ref = na.omit(est01.bs[, v]) - est01$zero[v]
    # summand = mean(abs(ref) >= abs(est01$zero[v]))
    ltpv = mean(ref <= est01$zero[v])
    summand = 2 * min(ltpv, 1 - ltpv)
    p01.bs = c(p01.bs, summand)
  }
  
  for (i in 1:16){
    out.table[[cohort]][3*(i-1) + 1, 1] = est01$zero[i]
    out.table[[cohort]][3*(i-1) + 2, 1] = est01$zero[i + 16]
    out.table[[cohort]][3*(i-1) + 3, 1] = est01$zero[i + 2*16]
    
    out.table[[cohort]][3*(i-1) + 1, 2] = sqrt(-diag(solve(est01$hess)))[i]
    out.table[[cohort]][3*(i-1) + 2, 2] = sqrt(-diag(solve(est01$hess)))[i+16]
    out.table[[cohort]][3*(i-1) + 3, 2] = sqrt(-diag(solve(est01$hess)))[i+32]
    
    out.table[[cohort]][3*(i-1) + 1, 3] = se.bs[i]
    out.table[[cohort]][3*(i-1) + 2, 3] = se.bs[i+16]
    out.table[[cohort]][3*(i-1) + 3, 3] = se.bs[i+32]
    
    out.table[[cohort]][3*(i-1) + 1, 4] = 2 * (1 - pnorm(abs(out.table[[cohort]][3*(i-1) + 1, 1]/out.table[[cohort]][3*(i-1) + 1, 2]), 0, 1)) 
    out.table[[cohort]][3*(i-1) + 2, 4] = 2 * (1 - pnorm(abs(out.table[[cohort]][3*(i-1) + 2, 1]/out.table[[cohort]][3*(i-1) + 2, 2]), 0, 1))
    out.table[[cohort]][3*(i-1) + 3, 4] = 2 * (1 - pnorm(abs(out.table[[cohort]][3*(i-1) + 3, 1]/out.table[[cohort]][3*(i-1) + 3, 2]), 0, 1)) 
    
    out.table[[cohort]][3*(i-1) + 1, 5] = p01.bs[i]
    out.table[[cohort]][3*(i-1) + 2, 5] = p01.bs[i+16]
    out.table[[cohort]][3*(i-1) + 3, 5] = p01.bs[i+32]
  }
  
  out.table[[cohort]][49:(48+G), 1] = est01$zero[49:(48+G)]
  out.table[[cohort]][49:(48+G), 2] = sqrt(-diag(solve(est01$hess)))[49:(48+G)]
  out.table[[cohort]][49:(48+G), 3] = se.bs[49:(48+G)]
  out.table[[cohort]][49:(48+G), 4] = 2 * (1 - pnorm(abs(out.table[[cohort]][49:(48+G), 1]/out.table[[cohort]][49:(48+G), 2]), 0, 1)) 
  out.table[[cohort]][49:(48+G), 5] = p01.bs[49:(48+G)]
  
  order01 = match(est01$zero[p01.bs <= 0.05], est01$zero)
  order01.1 = order01[order01<=k]
  order01.2 = order01[(order01<=(2*k)) & (order01>k)]
  order01.3 = order01[(order01>(2*k + 1)) & (order01 <= (3*k))]
  coef01.1 = est01$zero[order01.1]
  coef01.2 = est01$zero[order01.2]
  coef01.3 = est01$zero[order01.3]
  
  x01.1 = X[[cohort]][, order01.1]
  x01.2 = X[[cohort]][, order01.2 - k]
  x01.3 = X[[cohort]][, order01.3 - 2*k]
  a01 = est01$zero[(3*k+1):(3*k+G)]
  
  P01 = matrix(0, nc=4, nr = n[cohort])
  for (i in 1:n[cohort]){
    
    a = ifelse( p01.bs[membership[i] + 3*k] > 0.05, 0, a01[membership[i]])
    
    eta1 = as.numeric(x01.1[i,] %*% coef01.1 + a)
    eta2 = as.numeric(x01.2[i,] %*% coef01.2 + a)
    eta3 = as.numeric(x01.3[i,] %*% coef01.3 + a)
    
    pi = c((1 - plogis(eta1)), (plogis(eta1)*(1 - plogis(eta2))), (plogis(eta1) * plogis(eta2) * (1 - plogis(eta3))), (plogis(eta1) * plogis(eta2) * plogis(eta3)))
    
    P01[i,] = pi
    
  }
  
  pred01 = rep(0, n[cohort]); correct01 = pred01; actu01 = rep(0, n[cohort])
  indicator = as.data.frame(cbind(dummy[[cohort]][,1], dummy[[cohort]][,2], dummy[[cohort]][,3], dummy[[cohort]][,4]))
  
  for (i in 1:n[cohort]){
    pred01[i] = match(max(P01[i,]), P01[i,])
    actu01[i] = match(max(indicator[i,]), indicator[i,])
    correct01[i] = ifelse(pred01[i] == actu01[i], 1, 0)
  }
  
  precision = sum(correct01==1)/length(correct01)
  accuracy[cohort, 2] = precision
  accuracy[cohort, 3] = n.equiv
  
  tab01 = cbind(actu01, pred01)
  confusion01 = matrix(0, nc=4, nr=4)
  for (j in 1:4){
    for (l in 1:4){
      njk = 0
      for (i in 1:n[cohort]){
        njk = njk + ifelse(((tab01[i,1]==j) & (tab01[i,2]==l)), 1, 0)
      }
      confusion01[j,l] = njk
    }
  }
  confusion[[cohort]] = confusion01 / apply(confusion01, 1, sum)
  
}

for (cohort in 1:2){
  write.csv(out.table[[cohort]], paste("model output cohort", cohort, ".csv"))
  write.csv(confusion[[cohort]], paste("confusion matrix cohort", cohort, ".csv"))
}


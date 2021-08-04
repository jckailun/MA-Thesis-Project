B = 1000
set.seed(1000)

n = c(nrow(X1.c), nrow(X2.c))
G = 1
H = c(2, 7)
LR.comp = rep(0, 2)
LR.bs.comp = rep(list(NULL), 2)
LR.dist.comp = rep(list(NULL), 2)
pval.lr.comp = rep(0, 2)

for (cohort in 1:2){
  
  theta0 = rep(0, 3*k + G)
  theta1 = rep(0, 3*k + H[cohort])
  
  X = list(X1.c, X2.c)
  X.c = list(X1.cc, X2.cc)
  diss = list(diss1, diss2)
  dummy = list(dummy.1, dummy.2)
  
  dat = clust(cohort, H[cohort], BS = FALSE)
  df = dat[[1]]
  D = dat[[2]]
  membership = dat[[3]]
  LR1 = brd(theta1, H[cohort], k=16, cohort = cohort)
  
  dat = clust(cohort, G, BS = FALSE)
  df = dat[[1]]
  D = dat[[2]]
  membership = dat[[3]]
  LR0 = brd(theta0, G, k=16, cohort = cohort)
  
  LR.comp[cohort] = - 2 * (LR0$loglik[length(LR0$loglik)] - LR1$loglik[length(LR1$loglik)])

  for (b in 1:B){
    
    idx.bs = sample(1:n[cohort], n[cohort], replace = TRUE)
    
    dummy.bs = dummy
    dummy.bs[[cohort]] = dummy.bs[[cohort]][idx.bs, ]
    
    theta.bs = rep(0, 3*k + H[cohort])
    a.bs = rep(0, 3*k + G)
    
    X.bs = X
    X.bs[[cohort]] = X.bs[[cohort]][idx.bs, ]
    
    X.bs.c = X.c
    X.bs.c[[cohort]] = X.bs.c[[cohort]][idx.bs, ]
    
    gower_dist <- daisy(X.bs.c[[cohort]], metric = "gower", type = list(symm = c(12, 14:17)))
    diss.bs = as.matrix(gower_dist)
    
    dat = clust(cohort, H[cohort], BS = TRUE)
    df = dat[[1]]
    D = dat[[2]]
    # membership = dat[[3]]
    LR1.bs = brd(theta.bs, H[cohort], k = 16, cohort = cohort)
    
    dat = clust(cohort, G, BS = TRUE)
    df = dat[[1]]
    D = dat[[2]]
    # membership = list(dat[[3]])
    LR0.bs = brd(a.bs, G, k = 16, cohort = cohort)
    
    LR.bs.comp[[cohort]] = rbind(LR.bs.comp[[cohort]], 
                            (- 2 * (LR0.bs$loglik[length(LR0.bs$loglik)] - LR1.bs$loglik[length(LR1.bs$loglik)]))
    )
    print((- 2 * (LR0.bs$loglik[length(LR0.bs$loglik)] - LR1.bs$loglik[length(LR1.bs$loglik)])))
  }
  
  LR.dist.comp[[cohort]] = LR.bs.comp[[cohort]] - LR.comp[cohort]
  pval.lr.comp[cohort] = mean(LR.dist.comp[[cohort]] >= LR.comp[cohort])
  
}

for (cohort in 1:2){
  A = LR.dist.comp[[cohort]][(length(LR.dist.comp[[cohort]])-999):length(LR.dist.comp[[cohort]])]
  B = LR.comp[cohort]
  plot(density(A)); abline(v = B)
  pval.lr.comp[cohort] = mean(A >= B)
}




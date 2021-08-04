B = 1000
set.seed(1000)
require(pracma)

n = c(nrow(X1.c), nrow(X2.c))
H = c(2, 7)
LR = rep(0, 2); mcfadden = LR
LR.bs = rep(list(NULL), 2)
LR.dist = rep(list(NULL), 2)
pval.lr = rep(0, 2)

for (cohort in 1:2){
  
  X = list(X1.c, X2.c)
  X.c = list(X1.cc, X2.cc)
  diss = list(diss1, diss2)
  dummy = list(dummy.1, dummy.2)
  
  G = H[cohort]
  theta0 = rep(0, 3*k + G)
  a0 = rep(0, G)
  
  dat = clust(cohort, G, BS = FALSE)
  df = dat[[1]]
  D = dat[[2]]
  membership = dat[[3]]
  
  LR1 = brd(theta0, G, k=16, cohort = cohort)
  LR0 = brd.lr(a0, cohort = cohort)
  LR[cohort] = - 2 * (LR0$loglik[length(LR0$loglik)] - LR1$loglik[length(LR1$loglik)])
  mcfadden[cohort] = 1 - (LR1$loglik[length(LR1$loglik)] / LR0$loglik[length(LR0$loglik)] )
  
  for (b in 1:B){

    idx.bs = sample(1:n[cohort], n[cohort], replace = TRUE)

    dummy.bs = dummy
    dummy.bs[[cohort]] = dummy.bs[[cohort]][idx.bs, ]

    theta.bs = rep(0, 3*k + G)
    a.bs = rep(0, G)
    
    X.bs = X
    X.bs[[cohort]] = X.bs[[cohort]][idx.bs, ]

    X.bs.c = X.c
    X.bs.c[[cohort]] = X.bs.c[[cohort]][idx.bs, ]

    gower_dist <- daisy(X.bs.c[[cohort]], metric = "gower", type = list(symm = c(12, 14:17)))
    diss.bs = as.matrix(gower_dist)

    dat = clust(cohort, G, BS = TRUE)

    df = dat[[1]]
    D = dat[[2]]
    
    LR1.bs = brd(theta.bs, G, k = 16, cohort = cohort)
    LR0.bs = brd.lr(a.bs, cohort = cohort)
    
    LR.bs[[cohort]] = rbind(LR.bs[[cohort]], 
    (- 2 * (LR0.bs$loglik[length(LR0.bs$loglik)] - LR1.bs$loglik[length(LR1.bs$loglik)]))
    )
    
  }
  
  LR.dist[[cohort]] = LR.bs[[cohort]] - LR[cohort]
  pval.lr[cohort] = mean(LR.dist[[cohort]] >= LR[cohort])
    
}
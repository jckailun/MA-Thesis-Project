B = 1000
set.seed(1000)

H = c(2,7)
X = list(X1.c, X2.c)
diss = list(diss1, diss2)
dummy = list(dummy.1, dummy.2)
ME.class = rep(list(matrix(0, nc = 12, nr = 16)), 2)
pval.me1 = rep(list(matrix(0, nr=16, nc=4)), 2)

for (cohort in 1:2){
  
  G = H[cohort]
  
  X = list(X1.c, X2.c); X.c = list(X1.cc, X2.cc)
  diss = list(diss1, diss2)
  dummy = list(dummy.1, dummy.2)
  
  dat = clust(cohort, G, BS = FALSE)
  df = dat[[1]]
  D = dat[[2]]
  
  theta0 = rep(0, 48 + G)
  
  est.me = brd(theta0, G, k = 16, cohort = cohort)
  b1 = list(est.me$zero[1:16], est.me$zero[17:32], est.me$zero[33:48])
  a = out$zero[49:(48+G)]
  
  ### ME calculation
  membership = dat[[3]]
  a1 = rep(0, n[cohort])
  for (i in 1:n[cohort]){ 
    a1[i] = a[membership[i]] 
  }
  for (k in 1:16){
    if (k %in% c(11, 13:16)){ 
      
      dat.temp = X[[cohort]]
      
      dat.temp[, k] = 1
      eta1 = dat.temp %*% b1[[1]] + a1
      eta12 = dat.temp %*% b1[[2]] + a1
      eta13 = dat.temp %*% b1[[3]] + a1
      
      dat.temp[, k] = 0
      eta0 = dat.temp %*% b1[[1]] + a1
      eta02 = dat.temp %*% b1[[2]] + a1
      eta03 = dat.temp %*% b1[[3]] + a1
      
      ME.class[[cohort]][k, 1] = mean((1-plogis(eta1)) - (1-plogis(eta0)))
      
      ME.class[[cohort]][k, 4] = mean(plogis(eta1)*(1-plogis(eta12)) - plogis(eta0)*(1-plogis(eta02)))
      
      ME.class[[cohort]][k, 7] = mean(plogis(eta1)*(plogis(eta12)*(1-plogis(eta13))) - plogis(eta0)*plogis(eta02)*(1-plogis(eta03)))
      
      ME.class[[cohort]][k, 10] = mean(plogis(eta1)*(plogis(eta12)*plogis(eta13)) - plogis(eta0)*plogis(eta02)*plogis(eta03))
      
    } else {
      
      eta1 = X[[cohort]] %*% b1[[1]] + a1
      eta2 = X[[cohort]] %*% b1[[2]] + a1
      eta3 = X[[cohort]] %*% b1[[3]] + a1
      
      ME.class[[cohort]][k, 1] = mean(  - plogis(eta1) * (1-plogis(eta1)) * b1[[1]][k] )
      ME.class[[cohort]][k, 4] = mean(  (1 - plogis(eta2)) * plogis(eta1) * (1-plogis(eta1)) *  b1[[1]][k] +
                                 plogis(eta1) * (- plogis(eta2) * (1-plogis(eta2))) * b1[[2]][k]  ) 
      ME.class[[cohort]][k, 7] = mean( (1 - plogis(eta3)) * plogis(eta2) * (plogis(eta1) * (1-plogis(eta1))) * b1[[1]][k] +
                                (1 - plogis(eta3)) * plogis(eta1) * (plogis(eta2) * (1-plogis(eta2))) * b1[[2]][k] + 
                                plogis(eta1) * plogis(eta2) * (-plogis(eta3) * (1-plogis(eta3))) * b1[[3]][k] )
      ME.class[[cohort]][k, 10] = mean( plogis(eta3) * plogis(eta2) * plogis(eta1) * (1-plogis(eta1)) * b1[[1]][k] +
                                 plogis(eta3) * plogis(eta1) * plogis(eta2) * (1-plogis(eta2)) * b1[[2]][k] + 
                                 plogis(eta1) * plogis(eta2) * plogis(eta3) * (1-plogis(eta3)) * b1[[3]][k] ) 
      
    }
  }
  ### End of ME calculation ###

  # idj = list(1:3, 4:6, 7:9, 10:12)
  
  ### ME bootstrapping ###
  
  ME.bs.c = list(matrix(0, nc=16, nr=B), matrix(0, nc=16, nr=B), matrix(0, nc=16, nr=B), matrix(0, nc=16, nr=B))
  
  for (b in 1:B){
  
    idx.bs = sample(1:n[cohort], n[cohort], replace = TRUE)
    dummy.bs = dummy
    dummy.bs[[cohort]] = dummy.bs[[cohort]][idx.bs, ]
  
    theta.bs = rep(0, 3*k + G)
  
    X.bs = X
    X.bs[[cohort]] = X.bs[[cohort]][idx.bs, ]
  
    X.bs.c = X.c
    X.bs.c[[cohort]] = X.bs.c[[cohort]][idx.bs, ]
  
    # reshuffle the data set 
  
    gower_dist <- daisy(X.bs.c[[cohort]], metric = "gower", type = list(symm = c(12, 14:17)))
    diss.bs = as.matrix(gower_dist)
    dat = clust(cohort, G, BS = TRUE)
    # setting BS = T makes dat absorb and utilise .bs objects
    # and then cluster
  
    df = dat[[1]]
    D = dat[[2]]
    # since df and D are entirely controlled by G and cohort,
    # it does not matter if we overwrote them in the BS procedure.
    # en fin de la journee, at the beginning of each iteration, 
    # we re-adjust G. For each G, we have one BS run. 
    membership = dat[[3]]
  
    out = brd(theta.bs, G, k = 16, cohort = cohort)
    # est01.bs = rbind(est01.bs, out$zero)
    # finally obtain the estimates using the clustered data set
  
    b1.bs = list(out$zero[1:16], out$zero[17:32], out$zero[33:48])
    a1.bs = out$zero[49:(48+G)]
    a.bs = rep(0, n[cohort])
    for (i in 1:n[cohort]){ 
      a.bs[i] = a1.bs[membership[i]] 
    }
    
    for (k in 1:16){
      if (k %in% c(11, 13:16)){
      
        dat.temp = X.bs[[cohort]]
      
        dat.temp[, k] = 1
        eta1.bs = dat.temp %*% b1.bs[[1]] + a.bs
        eta12.bs = dat.temp %*% b1.bs[[2]] + a.bs
        eta13.bs = dat.temp %*% b1.bs[[3]] + a.bs
      
        dat.temp[, k] = 0
        eta0.bs = dat.temp %*% b1.bs[[1]] + a.bs
        eta02.bs = dat.temp %*% b1.bs[[2]] + a.bs
        eta03.bs = dat.temp %*% b1.bs[[3]] + a.bs
      
        ME.bs.c[[1]][b, k] = mean((1-plogis(eta1.bs)) - (1-plogis(eta0.bs)))
      
        ME.bs.c[[2]][b, k] = mean(plogis(eta1.bs)*(1-plogis(eta12.bs)) - plogis(eta0.bs)*(1-plogis(eta02.bs)))
      
        ME.bs.c[[3]][b, k] = mean(plogis(eta1.bs)*(plogis(eta12.bs)*(1-plogis(eta13.bs))) - plogis(eta0.bs)*plogis(eta02.bs)*(1-plogis(eta03.bs)))
      
        ME.bs.c[[4]][b, k] = mean(plogis(eta1.bs)*(plogis(eta12.bs)*plogis(eta13.bs)) - plogis(eta0.bs)*plogis(eta02.bs)*plogis(eta03.bs))
      
      } else {
      
        eta1.bs = X.bs[[cohort]] %*% b1.bs[[1]] + a.bs
        eta2.bs = X.bs[[cohort]] %*% b1.bs[[2]] + a.bs
        eta3.bs = X.bs[[cohort]] %*% b1.bs[[3]] + a.bs
      
        ME.bs.c[[1]][b, k] = mean(  - plogis(eta1.bs) * (1-plogis(eta1.bs)) * b1.bs[[1]][k] )
        ME.bs.c[[2]][b, k] = mean(  (1 - plogis(eta2.bs)) * plogis(eta1.bs) * (1-plogis(eta1.bs)) *  b1.bs[[1]][k] +
                                    plogis(eta1.bs) * (- plogis(eta2.bs) * (1-plogis(eta2.bs))) * b1.bs[[2]][k]  ) 
        ME.bs.c[[3]][b, k] = mean( (1 - plogis(eta3.bs)) * plogis(eta2.bs) * (plogis(eta1.bs) * (1-plogis(eta1.bs))) * b1.bs[[1]][k] +
                                   (1 - plogis(eta3.bs)) * plogis(eta1.bs) * (plogis(eta2.bs) * (1-plogis(eta2.bs))) * b1.bs[[2]][k] + 
                                   plogis(eta1.bs) * plogis(eta2.bs) * (-plogis(eta3.bs) * (1-plogis(eta3.bs))) * b1.bs[[3]][k] )
        ME.bs.c[[4]][b, k] = mean( plogis(eta3.bs) * plogis(eta2.bs) * plogis(eta1.bs) * (1-plogis(eta1.bs)) * b1.bs[[1]][k] +
                                   plogis(eta3.bs) * plogis(eta1.bs) * plogis(eta2.bs) * (1-plogis(eta2.bs)) * b1.bs[[2]][k] + 
                                   plogis(eta1.bs) * plogis(eta2.bs) * plogis(eta3.bs) * (1-plogis(eta3.bs)) * b1.bs[[3]][k] ) 
      
      }
    }
  }
  
  ### End of ME bootstrapping ###
  
  se.c1 = NULL
  for (j in 1:4){
    se.c1 = cbind(se.c1, sqrt(diag(cov(ME.bs.c[[j]]))))
  }

  ME.class[[cohort]][,2] = se.c1[,1]; ME.class[[cohort]][,3] = ME.class[[cohort]][,1] / ME.class[[cohort]][,2]
  ME.class[[cohort]][,5] = se.c1[,2]; ME.class[[cohort]][,6] = ME.class[[cohort]][,4] / ME.class[[cohort]][,5]
  ME.class[[cohort]][,8] = se.c1[,3]; ME.class[[cohort]][,9] = ME.class[[cohort]][,7] / ME.class[[cohort]][,8]
  ME.class[[cohort]][,11] = se.c1[,4]; ME.class[[cohort]][,12] = ME.class[[cohort]][,10] / ME.class[[cohort]][,11]

  ref = ME.class[[cohort]][,c(1,4,7,10)]
  for (j in 1:4){
    for (k in 1:16){
       ltpv = mean(abs(ME.bs.c[[j]][,k]-ref[,j][k]) <= abs(ref[,j][k]))
       summand = 2 * min(ltpv, 1 - ltpv)
       pval.me1[[cohort]][k, j] = summand
    }
  }

}

ME.output = rep(list(matrix(0, nr=16, nc=16)), 2)
for (cohort in 1:2){
  ME.output[[cohort]][,1:3] = ME.class[[cohort]][,1:3]
  ME.output[[cohort]][,5:7] = ME.class[[cohort]][,4:6]
  ME.output[[cohort]][,9:11] = ME.class[[cohort]][,7:9]
  ME.output[[cohort]][,13:15] = ME.class[[cohort]][,10:12]
  ME.output[[cohort]][,c(4,8,12,16)] = pval.me1[[cohort]]
  ME.output[[cohort]] = as.data.frame(ME.output[[cohort]])
  names(ME.output[[cohort]]) = rep(c('ME', 'Standard Error', 'Z', 'p-value'), 4)
  #write.csv(ME.output[[cohort]], paste("Marginal Effects of 4 Levels with BS cohort", cohort, ".csv"))
  write.table(ME.output[[cohort]], paste("Marginal Effects of 4 Levels with BS cohort", cohort, ".xls"))
} 

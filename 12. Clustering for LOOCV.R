clust.loocv = function(cohort, G, BS = FALSE){

  if (BS == FALSE){
    dis = diss.loocv
    dat = X.seen
    dum = dummy.seen
  } 
  if (BS == TRUE) { 
    dis = diss.loocv.bs
    dat = X.loocv.bs
    dum = dummy.loocv.bs
  }

  n = nrow(dis)
    
  if (G == 1){
    grp = list(dat)
    d = list(dum)
    clus = rep(1, n)
    D = NULL
    E = NULL
  } else {
      
  aaa = pam(dis, G, diss = TRUE)
  clus = aaa$clustering
  A = aaa$clusinfo
  B = aaa$silinfo$clus.avg.widths 
  C = cbind(A, B)
  D = aaa$medoids
  E = aaa$id.med
      
  # returns the estimated cluster membership of each observation
  idx = rep(list(0), G)
  grp = rep(list(0), G)
  d = rep(list(0), G)
      
  for (j in 1:G){
        
    for (i in 1:n){
      idx[[j]] = c(idx[[j]], ifelse(clus[i] == j, i, 0))
    }
      idx[[j]] = idx[[j]][idx[[j]] != 0]
        
      grp[[j]] = dat[idx[[j]], ]
      d[[j]] = dum[idx[[j]], ]
        
   }
  }
    
  list(grpeddata = grp, dummy = d, membership = clus, clusinfo = C, med = D, med2 = E) 
    
}


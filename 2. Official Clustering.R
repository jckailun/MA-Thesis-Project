# Clustering with PAM #
install.packages("ClusterR"); library(ClusterR)
install.packages("cluster"); library(cluster)

gower_dist_1 <- daisy(X1.cc, metric = "gower", type = list(symm = c(12, 14:17)))
gower_dist_2 <- daisy(X2.cc, metric = "gower", type = list(symm = c(12, 14:17)))
diss1 = as.matrix(gower_dist_1)
diss2 = as.matrix(gower_dist_2)

X = list(X1.c, X2.c)
diss = list(diss1, diss2)

dummy.1 = cbind(D11, D21, D31, D41)
dummy.2 = cbind(D12, D22, D32, D42)
dummy = list(dummy.1, dummy.2)

clust = function(cohort, G, BS = FALSE){
  
  if (BS == FALSE){
    dis = diss[[cohort]]
    dat = X[[cohort]]
    dum = dummy[[cohort]]
    } 
  if (BS == TRUE) { 
    dis = diss.bs
    # the quantity diss.bs comes from X.bs.c - which includes acgrd
    dat = X.bs[[cohort]]
    dum = dummy.bs[[cohort]]
  }

  n = nrow(dis)
  
  if (G == 1){
    grp = list(dat)
    d = list(dum)
    clus = rep(1, n)
    } else {
  
  aaa = pam(dis, G, diss = TRUE)
  clus = aaa$clustering
  A = aaa$clusinfo
  B = aaa$silinfo$clus.avg.widths
  C = cbind(A, B)
  D = aaa$medoids
  
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
  
  list(grp, d, clus, clusinfo = C, D) 
  
}

X = list(X1.c, X2.c)
diss = list(diss1, diss2)
dummy = list(dummy.1, dummy.2)
H = c(3, 7)
for (cohort in 1:2){
  
  G = H[cohort]
  info = clust(cohort, G, BS = FALSE)$clusinfo
  info = as.data.frame(info)
  write.csv(info, paste("Cluster Info cohort", cohort, ".csv"))

}
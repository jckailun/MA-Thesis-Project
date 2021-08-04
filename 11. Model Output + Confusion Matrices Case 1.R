##### Output Table + Confusion Matrix w/ G = 1 #####
require(MASS)

set.seed(1000)

k = 16
H = rep(1, 2)
n = c(nrow(X1.c), nrow(X2.c))
B = 1000

accuracy = as.data.frame(cbind(H, rep(0, length(H)), rep(0, length(H))))
names(accuracy) = c("G", "Model Precision", "Equivalent Bootstrap Size")
out.table.2 = list(matrix(0, nr = 3*k+H[1], nc = 3), matrix(0, nr = 3*k+H[2], nc = 3))
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
  se01 = sqrt(diag(-solve(est01$hess)))
  pv01 = 2 * (1 - pnorm(abs(est01$zero/se01), 0, 1))
    
  for (i in 1:16){
    out.table.2[[cohort]][3*(i-1) + 1, 1] = est01$zero[i]
    out.table.2[[cohort]][3*(i-1) + 2, 1] = est01$zero[i + 16]
    out.table.2[[cohort]][3*(i-1) + 3, 1] = est01$zero[i + 2*16]
    
    out.table.2[[cohort]][3*(i-1) + 1, 2] = sqrt(-diag(solve(est01$hess)))[i]
    out.table.2[[cohort]][3*(i-1) + 2, 2] = sqrt(-diag(solve(est01$hess)))[i+16]
    out.table.2[[cohort]][3*(i-1) + 3, 2] = sqrt(-diag(solve(est01$hess)))[i+32]
    
    out.table.2[[cohort]][3*(i-1) + 1, 3] = 2 * (1 - pnorm(abs(out.table.2[[cohort]][3*(i-1) + 1, 1]/out.table.2[[cohort]][3*(i-1) + 1, 2]), 0, 1)) 
    out.table.2[[cohort]][3*(i-1) + 2, 3] = 2 * (1 - pnorm(abs(out.table.2[[cohort]][3*(i-1) + 2, 1]/out.table.2[[cohort]][3*(i-1) + 2, 2]), 0, 1))
    out.table.2[[cohort]][3*(i-1) + 3, 3] = 2 * (1 - pnorm(abs(out.table.2[[cohort]][3*(i-1) + 3, 1]/out.table.2[[cohort]][3*(i-1) + 3, 2]), 0, 1)) 
    
  }
  out.table.2[[cohort]][49:(48+G), 1] = est01$zero[49:(48+G)]
  out.table.2[[cohort]][49:(48+G), 2] = sqrt(-diag(solve(est01$hess)))[49:(48+G)]
  out.table.2[[cohort]][49:(48+G), 3] = 2 * (1 - pnorm(abs(out.table.2[[cohort]][49:(48+G), 1]/out.table.2[[cohort]][49:(48+G), 2]), 0, 1)) 

  order01 = match(est01$zero[pv01 <= 0.05], est01$zero)
  # order01 = match(est01$zero[out.table.2[[cohort]][,3] <= 0.05], est01$zero) 
  order01.1 = order01[order01 <= k]
  order01.2 = order01[(order01 <= (2*k)) & (order01 > k)]
  order01.3 = order01[(order01 > (2*k + 1)) & (order01 <= (3*k))]
  coef01.1 = est01$zero[order01.1]
  coef01.2 = est01$zero[order01.2]
  coef01.3 = est01$zero[order01.3]
  
  x01.1 = X[[cohort]][, order01.1]
  x01.2 = X[[cohort]][, order01.2 - k]
  x01.3 = X[[cohort]][, order01.3 - 2*k]
  a01 = est01$zero[(3*k+1):(3*k+G)]
  
  P01 = matrix(0, nc=4, nr = n[cohort])
  for (i in 1:n[cohort]){
    
    a = ifelse( out.table.2[[cohort]][,3][length(theta0)] > 0.05, 0, a01)
    
    eta1 = as.numeric(x01.1[i,] %*% coef01.1 + a)
    eta2 = as.numeric(x01.2[i,] %*% coef01.2 + a)
    eta3 = as.numeric(x01.3[i,] %*% coef01.3 + a)
    
    pi = c((1 - plogis(eta1)), (plogis(eta1)*(1 - plogis(eta2))), (plogis(eta1) * plogis(eta2) * (1 - plogis(eta3))), (plogis(eta1) * plogis(eta2) * plogis(eta3)))
    
    P01[i,] = pi
    
  }
  
  pred01 = rep(0, n[cohort]); correct01 = rep(0, n[cohort]); actu01 = rep(0, n[cohort])
  indicator = as.data.frame(cbind(dummy[[cohort]][,1], dummy[[cohort]][,2], dummy[[cohort]][,3], dummy[[cohort]][,4]))
  
  for (i in 1:n[cohort]){
    pred01[i] = match(max(P01[i,]), P01[i,])
    actu01[i] = match(max(indicator[i,]), indicator[i,])
    correct01[i] = ifelse(pred01[i] == actu01[i], 1, 0)
  }
  
  precision = sum(correct01==1)/length(correct01)
  accuracy[cohort, 2] = precision

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
  write.csv(out.table.2[[cohort]], paste("model output2 cohort", cohort, " case 1.csv"))
  write.csv(confusion[[cohort]], paste("confusion matrix cohort", cohort, "case 1.csv"))
}


### Model Precision Output ###
precision.g_1 = accuracy$`Model Precision`*100
cohort1 = as.data.frame(cbind(1:20, c(precision.g_1[1],
49.4433,
53.3907,
52.2267,
52.7328,
51.1640,
47.9251,
46.5587,
46.7611,
46.6599,
47.5202,
50.3036,
46.5587,
45.1417,
47.0648,
43.7247,
43.5223,
46.4575,
43.4717,
42.4089)))
names(cohort1) = c("G", "Model Precision")

cohort2 = as.data.frame(cbind(1:20, c(precision.g_1[2], 
43.5683,
44.0462,
36.9574,
39.4265,
36.6388,
45.0418,
42.6125,
40.2628,
40.1832,
40.8602,
32.4174,
30.7447,
34.1298,
37.4353,
37.1565,
24.5321,
24.0542,
24.2135,
24.2533)))
names(cohort2) = c("G", "Model Precision")

par(mfrow = c(2,1))

plot(cohort1$G, cohort1$`Model Precision`, xlab = "G", ylab = "Prediction Accuracy (%)", main = "Cohort 1", ylim = c(0, 80))
abline(v = match(max(cohort1$`Model Precision`), cohort1$`Model Precision`), lty=3, lwd=2, col="Blue")
plot(cohort2$G, cohort2$`Model Precision`, xlab = "G", ylab = "Prediction Accuracy (%)", main = "Cohort 2", ylim = c(0, 80))
abline(v = match(max(cohort2$`Model Precision`), cohort2$`Model Precision`), lty=3, lwd=2, col="Blue")

#write.csv(cohort1, "model selection cohort 1.csv")
#write.csv(cohort2, "model selection cohort 2.csv")

for (cohort in 1:2){
  write.csv(confusion[[cohort]], paste("confusion matrix cohort", cohort, ".csv"))
}
  
  
  
### Output Table ###
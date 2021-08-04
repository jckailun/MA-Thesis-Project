##### Estimation with Hyperparameter G > 1 : Case 2 #####
require(MASS)
require(cluster)
require(ClusterR)

set.seed(1000)

B = 1050
H = 1:20
prec = rep(list(as.data.frame(cbind(H, rep(0, length(H)), rep(0, length(H)), 
                                    rep(0, length(H)), rep(0, length(H)), rep(0, length(H)), rep(0, length(H)), rep(0, length(H)) ))), 2)
names(prec[[1]]) = c("G", "Model Precision", "Equivalent Bootstrap Size", "Brier Score", "Entropy_2", "Entropy_10", "Entropy_e")
names(prec[[2]]) = c("G", "Model Precision", "Equivalent Bootstrap Size", "Brier Score", "Entropy_2", "Entropy_10", "Entropy_e")

confusion.selection = list(rep(list(matrix(0, nc = 4, nr = 4)), length(H)), rep(list(matrix(0, nc = 4, nr = 4)), length(H)))



k = 16
ptm <- proc.time()

for (cohort in 1:2){
  
  state_dist = apply(dummy[[cohort]], 2, mean)
  
for (z in 1:length(H)){

G = H[z]

X = list(X1.c, X2.c)
X.c = list(X1.cc, X2.cc)
n = c(nrow(X1.c), nrow(X2.c))
diss = list(diss1, diss2)
dummy = list(dummy.1, dummy.2)

dat = clust(cohort, G, BS = FALSE)
df = dat[[1]]
D = dat[[2]]
# this is the partitioned data set for dummies - by the clusters

theta0 = rep(0, 3*k + G)
membership = dat[[3]]

est01 = brd(theta0, G, k=16, cohort = cohort)
# this is the estimated parameters
# next we ought to bootstrap the standard error

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
  
  out = brd(theta.bs, G, k = 16, cohort = cohort)
  est01.bs = rbind(est01.bs, out$zero)
  # finally obtain the estimates using the clustered data set
}

mid = apply(est01.bs, 2, median)
mn = apply(est01.bs, 2, mean)

for (i in 1:B){
  if (abs(mean(est01.bs[i, ])) > max(mn)) { est01.bs[i, ] = NA }
}

se.bs = sqrt(diag(cov(na.omit(est01.bs)[1:1000,])))
z01 = est01$zero / se.bs
n.equiv = nrow(na.omit(est01.bs))

p01.bs = NULL #; CI.bs = NULL
for (v in 1:length(theta0)){
  ref = na.omit(est01.bs[, v])[1:1000] - est01$zero[v]
  ltpv = mean(ref <= est01$zero[v])
  summand = 2 * min(ltpv, 1 - ltpv)
  # summand = mean(abs(ref) >= abs(est01$zero[v]))
  p01.bs = c(p01.bs, summand)  ###
  
  # a = quantile(ref, 0.025)
  # b = quantile(ref, 0.975)
  # CI.bs = rbind(CI.bs, c(est01$zero[v] - b, est01$zero[v] - a))
}

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

pred01 = rep(0, n[cohort]); correct01 = rep(0, n[cohort]); actu01 = rep(0, n[cohort])
indicator = as.data.frame(cbind(dummy[[cohort]][,1], dummy[[cohort]][,2], dummy[[cohort]][,3], dummy[[cohort]][,4]))

for (i in 1:n[cohort]){
  pred01[i] = match(max(P01[i,]), P01[i,])
  actu01[i] = match(max(indicator[i,]), indicator[i,])
  correct01[i] = ifelse(pred01[i] == actu01[i], 1, 0)
}

precision = sum(correct01==1)/length(correct01)
prec[[cohort]][z, 2] = precision
prec[[cohort]][z, 3] = n.equiv

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

confusion.selection[[z]] = confusion01 / apply(confusion01, 1, sum)

brier = sum((P01 - dummy[[cohort]])^2)/n[cohort]
prec[[cohort]][z, 4] = brier

cross_entropy = -sum(state_dist %*% log2(t(P01)))/n[cohort]
prec[[cohort]][z, 5] = cross_entropy

cross_entropy.2 = -sum(state_dist %*% log10(t(P01)))/n[cohort]
prec[[cohort]][z, 6] = cross_entropy.2

cross_entropy.3 = -sum(state_dist %*% log(t(P01)))/n[cohort]
prec[[cohort]][z, 7] = cross_entropy.3

Hp = -sum(state_dist * log2(state_dist))
KL = cross_entropy - Hp
prec[[cohort]][z, 8] = KL
}
}

for (cohort in 1:2){
  write.csv(prec[[cohort]], paste("Overall Precision Cohort", cohort, ".csv"))
}
##### END Case 2 ####         

par(mfrow = c(2,1))
cohort1 = prec[[1]]
cohort2 = prec[[2]]

plot(cohort1$G, cohort1$`Model Precision`*100, xlab = "G", ylab = "Prediction Accuracy (%)", main = "Cohort 1", ylim = c(0, 80))
abline(v = match(max(cohort1$`Model Precision`), cohort1$`Model Precision`), lty=3, lwd=2, col="Blue")
plot(cohort2$G, cohort2$`Model Precision`*100, xlab = "G", ylab = "Prediction Accuracy (%)", main = "Cohort 2", ylim = c(0, 80))
abline(v = match(max(cohort2$`Model Precision`), cohort2$`Model Precision`), lty=3, lwd=2, col="Blue")

KL = matrix(0, nc = 2, nr = 20)
for (cohort in 1:2){
  state_dist = apply(dummy[[cohort]], 2, mean)
  Hp = -sum(state_dist * log2(state_dist))
  KL[, cohort] = prec[[cohort]][, 5] - Hp
  prec[[cohort]][, 8] = KL[, cohort]
  
  write.csv(prec[[cohort]], paste("model selection cohort", cohort, ".csv"))
}






### Conclusion: G1 = 3; G2 = 7 ###






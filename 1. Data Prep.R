### 1. DATA ###
require(haven); library(haven)

setwd('/Users/caokailun/Desktop/redo/Data')
d1 = read.csv("NLSY79.csv")
d2 = read.csv("NLSY97.csv")

vars1 = c("ASVAB_gs", "ASVAB_ar" , "ASVAB_wk" , "ASVAB_pc" , "ASVAB_no" , "ASVAB_cs", "ASVAB_mk", "ASVAB_mc", "ASVAB_ei", 
          "med", paste('acgrd', 1:13, sep = '_'), 'male', 'income', 'nuclear', 'urban', 'race')
vars1 = match(vars1, names(d1))
d1 = d1[, vars1]
black1 = ifelse(d1$race==1, 1, 0); hisp1 = ifelse(d1$race==2, 1, 0); others1 = ifelse(d1$race == 3, 1, 0)
d1 = as.data.frame(cbind(d1, black1, hisp1))
d1 = d1[,-28]

rowstodelete1 = NULL
for (i in 1:nrow(d1)){
  if (sum(is.na(d1[i,])) != 0){
    rowstodelete1 = c(rowstodelete1, i)}
}

d1 = d1[-rowstodelete1,]


vars2 = c("ASVAB_GS", "ASVAB_AR" , "ASVAB_WK" , "ASVAB_PC" , "ASVAB_NO" , "ASVAB_CS", "ASVAB_MK", "ASVAB_MC", 'ASVAB_EI', 
          "med", paste('acgrd', 1:13, sep = '_'), 'sex', 'income', 'nuclear', 'urban', 'race')
vars2 = match(vars2, names(d2))
d2 = d2[, vars2]
d2$race = ifelse(d2$race==4, 3, d2$race); d2$sex = ifelse(d2$sex==1, 1, 0)
black2 = ifelse(d2$race==1, 1, 0); hisp2 = ifelse(d2$race==2, 1, 0); others2 = ifelse(d2$race == 3, 1, 0)
d2 = as.data.frame(cbind(d2, black2, hisp2))
d2 = d2[,-28]

rowstodelete2 = NULL
for (i in 1:nrow(d2)){
  if (sum(is.na(d2[i,])) != 0){
    rowstodelete2 = c(rowstodelete2, i)}
}

d2 = d2[-rowstodelete2,]

# Next, let's compute the rates of having 0s' in ACGRDs
zeros1 = NULL
for (j in 11:23){ zeros1 = c(zeros1, sum((d1[,j]==0))) }
rate1 = zeros1/nrow(d1); margin1 = diff(zeros1)/nrow(d1)
# perhaps we shall choose asgrd_9 as the optimal grade transition variable.

zeros2 = NULL
for (j in 11:23){ zeros2 = c(zeros2, sum((d2[,j]==0))) }
rate2 = zeros2/nrow(d2); margin2 = diff(zeros2)/nrow(d2)
# perhaps we shall choose asgrd_9 as the optimal grade transition variable.

rowstodelete1 = NULL
for (i in 1:nrow(d1)){
  if (d1$acgrd_9[i] == 0){
    rowstodelete1 = c(rowstodelete1, i)}
}

rowstodelete2 = NULL
for (i in 1:nrow(d2)){
  if (d2$acgrd_10[i] == 0){
    rowstodelete2 = c(rowstodelete2, i)}
}

d1 = d1[-rowstodelete1,]
d2 = d2[-rowstodelete2,]

rowstodelete = NULL
for (i in 1:nrow(d1)){
  if (d1$med[i] == 0){
    rowstodelete = c(rowstodelete, i)}
}
d1 = d1[-rowstodelete, ]

n1 = nrow(d1)
n2 = nrow(d2)

d1[,1:9] = scale(d1[,1:9])
d2[,1:9] = scale(d2[,1:9])

### End of 1. DATA ###

### 2. Model Building ###
 
D11 = ifelse(d1$acgrd_9 <= 11, 1, 0) # dummy variable indicating that they did not transition at all from 1
D21 = ifelse(d1$acgrd_9 == 12, 1, 0) # dummy variable indicating that they transitioned from 1 - 2
D31 = ifelse(d1$acgrd_9 >= 13 & d1$acgrd_9 <= 15, 1, 0) # dummy variable indicating that they transitioned from 2 - 3
D41 = 1 - D11 - D21 - D31 # dummy variable indicating that they transitioned from 3 - 4

D12 = ifelse(d2$acgrd_10 <= 11, 1, 0) # dummy variable indicating that they did not transition at all from 1
D22 = ifelse(d2$acgrd_10 == 12, 1, 0) # dummy variable indicating that they transitioned from 1 - 2
D32 = ifelse(d2$acgrd_10 >= 13 & d2$acgrd_10 <= 15, 1, 0) # dummy variable indicating that they transitioned from 2 - 3
D42 = 1 - D12 - D22 - D32 # dummy variable indicating that they transitioned from 3 - 4

X1 = as.matrix(d1); acgrd1 = X1[,c(11:23)]; X1.c = X1[,-c(11:23)]
X2 = as.matrix(d2); acgrd2 = X2[,c(11:23)]; X2.c = X2[,-c(11:23)]
X1.c[, c(1:10, 12)] = scale(X1.c[, c(1:10, 12)])
X2.c[, c(1:10, 12)] = scale(X2.c[, c(1:10, 12)])

X1.cc = X1[,-c(11:18, 20:23)]
X2.cc = X2[,-c(11:19, 21:23)]
X1.cc[,c(10,11,13)] = scale(X1.cc[,c(10,11,13)])
X2.cc[,c(10,11,13)] = scale(X2.cc[,c(10,11,13)])

### End of 2. Model Building ###
###Here is a simulation study for case 1 linear models

#ensure empty environment
rm(list = ls())

#library required
library(doParallel)
library(MASS)
library(glmnet)

#set up working directory
setwd("/mnt/h/UbuntuRv2/Gibbs-sampler-algorithm/block1")

set.seed(101)

#simulation set up
#200 train samples and 200 test samples with 1000 predictors
n <- 200
p <- 1000

#correlation matrix M with rho_ij = rho^|i-j|
M <- diag(1,p)
rho <- 0.7
for (i in 1:p)
{
  for (j in 1:i)
  {
    M[j,i] <- rho^{i-j}
    M[i,j] <- M[j,i]
  }
}

### 100 times model averaging results
beta <- (1:15)/5

v.lasso <- vector()

j <- 0
while(j < 50){

  #data matrix
  x <- mvrnorm(n,rep(0,p),M)
  colnames(x) <- paste("x",1:p, sep = "")

  y <- x[,1:15]%*%beta + sin(rnorm(n)*pi) + cos(rnorm(n)*pi) + rnorm(n)

  m.cv <- cv.glmnet(x,y,family = "gaussian")
  t <- rep(0,p)
  t[as.vector(coef(m.cv)[-1])!=0] <- 1
  v.lasso <- rbind(v.lasso, t)

  j <- j+1
  print(j)
}

save.image("simu_case1_lasso_07.RData")

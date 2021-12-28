###Here is a simulation study for case 1 linear models

#ensure empty environment
rm(list = ls())

#library required
library(doParallel)
library(MASS)
library(SIS)

#set up working directory
setwd("/mnt/h/UbuntuRv2/Gibbs-sampler-algorithm/block1")

set.seed(101)

#simulation set up
#200 train samples and 200 test samples with 1000 predictors
n <- 500
p <- 1000

#correlation matrix M with rho_ij = rho^|i-j|
M <- diag(1,p)
rho <- 0.5
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

v.vani.bic <- vector()
v.vani.ebic <- vector()

v.aggr.bic <- vector()
v.aggr.ebic <- vector()

j <- 0
while(j < 50){

  #data matrix
  x <- mvrnorm(n,rep(0,p),M)
  colnames(x) <- paste("x",1:p, sep = "")

  w <- x[,1:15]%*%beta + sin(rnorm(n)*pi) + cos(rnorm(n)*pi)
  q <- exp(w)/(1+exp(w))
  y <- rbinom(n,1,q)


  model = SIS(x, y, family = "binomial", penalty = "SCAD", tune = "bic",
              nsis = 100, varISIS = "vanilla", seed = 9,standardize = FALSE)
  t <- rep(0,p)
  t[model$ix0] <- 1
  v.vani.bic <- rbind(v.vani.bic, t)

  model = SIS(x, y, family = "binomial", penalty = "SCAD", tune = "ebic",
              nsis = 100, varISIS = "vanilla", seed = 9,standardize = FALSE)
  t <- rep(0,p)
  t[model$ix0] <- 1
  v.vani.ebic <- rbind(v.vani.ebic, t)

  model = SIS(x, y, family = "binomial", penalty = "SCAD", tune = "bic",
              nsis = 100, varISIS = "aggr", seed = 9,standardize = FALSE)
  t <- rep(0,p)
  t[model$ix0] <- 1
  v.aggr.bic <- rbind(v.aggr.bic, t)

  model = SIS(x, y, family = "binomial", penalty = "SCAD", tune = "ebic",
              nsis = 100, varISIS = "aggr", seed = 9,standardize = FALSE)
  t <- rep(0,p)
  t[model$ix0] <- 1
  v.aggr.ebic <- rbind(v.aggr.ebic, t)

  j <- j+1
  print(j)
}

save.image("simu_case1_isis_log_05.RData")

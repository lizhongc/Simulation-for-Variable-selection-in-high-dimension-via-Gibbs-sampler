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
n <- 200
p <- 1000

#correlation matrix
#correlation between any two predictors is 0.5 except x4,x5
rho  <- 0.5
M1   <- rho + (1-rho)*diag(p/2)

#correlation between x4 and others is 1/sqrt(0.5) except x5
M1[,4]  <- 1/sqrt(2)
M1[4,]  <- 1/sqrt(2)
M1[4,4] <- 1

M2 <- diag(1,p/2)
for (i in 1:p/2)
{
  for (j in 1:i)
  {
    M2[j,i] <- rho^{i-j}
    M2[i,j] <- M2[j,i]
  }
}

### 100 times model averaging results
beta <- rep(2,10)
beta[4] <- -6/sqrt(2)

v.vani.bic <- vector()
v.vani.ebic <- vector()

v.aggr.bic <- vector()
v.aggr.ebic <- vector()

j <- 0
while(j < 50){

  #data matrix
  x1 <- mvrnorm(n,rep(0,p/2),M1)
  x2 <- mvrnorm(n,rep(0,p/2),M2)
  x <- cbind(x1,x2)
  colnames(x) <- paste("x",1:p, sep = "")

  y <- x[,c(1:4,1:6+p/2)]%*%beta + sin(rnorm(n)*pi) + cos(rnorm(n)*pi) + rnorm(n)


  model = SIS(x, y, family = "gaussian", penalty = "SCAD", tune = "bic",
              nsis = 100, varISIS = "vanilla", seed = 9,standardize = FALSE)
  t <- rep(0,p)
  t[model$ix0] <- 1
  v.vani.bic <- rbind(v.vani.bic, t)

  model = SIS(x, y, family = "gaussian", penalty = "SCAD", tune = "ebic",
              nsis = 100, varISIS = "vanilla", seed = 9,standardize = FALSE)
  t <- rep(0,p)
  t[model$ix0] <- 1
  v.vani.ebic <- rbind(v.vani.ebic, t)

  model = SIS(x, y, family = "gaussian", penalty = "SCAD", tune = "bic",
              nsis = 100, varISIS = "aggr", seed = 9,standardize = FALSE)
  t <- rep(0,p)
  t[model$ix0] <- 1
  v.aggr.bic <- rbind(v.aggr.bic, t)

  model = SIS(x, y, family = "gaussian", penalty = "SCAD", tune = "ebic",
              nsis = 100, varISIS = "aggr", seed = 9,standardize = FALSE)
  t <- rep(0,p)
  t[model$ix0] <- 1
  v.aggr.ebic <- rbind(v.aggr.ebic, t)

  j <- j+1
  print(j)
}

save.image("simu_case1_isis_2.RData")

###Here is a simulation study for case 1 linear models

#ensure empty environment
rm(list = ls())

#library required
library(doParallel)
library(MASS)
library(IBGS)

#set up working directory
setwd("/nfs/ms_home/home/ad/student.unimelb.edu.au/lizhongc/myGit/iterated-block-gibbs/simu")

#parallel setting
registerDoParallel(30)
set.seed(101)

#simulation set up
#200 train samples and 200 test samples with 1000 predictors
n <- 200
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

v.prob <- vector()

j <- 0
while(j < 50){

  #data matrix
  x <- mvrnorm(n,rep(0,p),M)
  colnames(x) <- paste("x",1:p, sep = "")

  y <- x[,1:15]%*%beta + sin(rnorm(n)*pi) + cos(rnorm(n)*pi) + rnorm(n)

  m.block <- BlockGibbsSampler(y, x, info = "BIC", family = gaussian())
  v.prob <- rbind(v.prob, m.block$v.prob)

  j <- j+1
}

save.image("simu_case1_ibgs_bic_05.RData")

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
n <- 500
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

v.prob <- vector()

j <- 0
while(j < 50){

  #data matrix
  x1 <- mvrnorm(n,rep(0,p/2),M1)
  x2 <- mvrnorm(n,rep(0,p/2),M2)
  x <- cbind(x1,x2)
  colnames(x) <- paste("x",1:p, sep = "")

  w <- x[,c(1:4,1:6+p/2)]%*%beta + sin(rnorm(n)*pi) + cos(rnorm(n)*pi)
  q <- exp(w)/(1+exp(w))
  y <- rbinom(n,1,q)

  m.block <- BlockGibbsSampler(y, x, info = "exBIC", gamma = 0.5, family = "binomial")

  j <- j+1
}

save.image("simu_case2_ibgs_exbic_2.RData")

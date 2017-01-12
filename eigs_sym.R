library(microbenchmark)
library(mgcv)
library(RSpectra)
library(irlba)

n <- 1000
A <- matrix(rnorm(n*n),n,n)
m <- A + t(A)
microbenchmark(mgcv::slanczos(m,k=10,nt=1),
               RSpectra::eigs_sym(m,k=10),
               irlba(m,nv=10),
               times = 10)
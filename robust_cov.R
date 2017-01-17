# Draw the ellipse (what does gif mean in terms of ellipse scaling ?)
# Check visiondummmy.com

library(pcadapt)
library(ggplot2)
dt <- as.matrix(read.table(system.file("extdata", "geno3pops.pcadapt", package = "pcadapt")))
x <- pcadapt(dt, K = 2)
gt <- 1:150
xmat <- x$zscores

mrob = function(h){median(h, na.rm = TRUE)}
srob = function(h){mad(h, na.rm = TRUE)}

my_ogk = function(x, m = mrob, s = srob){
  n <- nrow(x)
  p <- ncol(x)
  D <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p){
    D[i, i] <- s(x[, i])
  }
  ymat <- matrix(0, nrow = n, ncol = p)
  for (i in 1:n){
    for (j in 1:p)
      ymat[i, j] = x[i, j] / D[j, j]
  }
  U <- matrix(0, nrow = p, ncol = p)
  for (j in 1:p){
    for (k in 1:p){
      U[j, k] = 1/(4 * (s(ymat[, j] + ymat[, k])^2 - s(ymat[, j] - ymat[, k])^2))
    }
  }
  E <- eigen(U, symmetric = TRUE)$vectors
  V <- ymat %*% E
  lambda <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p){
    lambda[i, i] <- srob(V[, i])^2
  }
  mvec <- matrix(0, nrow = p, ncol = 1)
  for (i in 1:p){
    mvec[i, 1] <- mrob(V[, i])
  }
  mu_vec <- E %*% mvec
  sigmat <- E %*% lambda %*% t(E)
  return(list(center = D %*% mu_vec, cov = D %*% sigmat %*% t(D)))
}

get.ellipse.coord = function(x, p, ci = 0.95, n.pts = 1000){
  xcenter <- x$center
  xcov <- x$cov
  ep <- eigen(xcov, symmetric = TRUE)
  s <- qchisq(ci, df = p)
  a <- 2 * sqrt(s * ep$values[1]) #90%
  b <- 2 * sqrt(s * ep$values[2])
  alpha <- atan(ep$vectors[2, 1] / ep$vectors[1, 1])
  t <- seq(0, 2 * pi, length.out = n.pts)
  x_coord <- a * cos(t) 
  y_coord <- b * sin(t) 
  coord <- as.matrix(rbind(x_coord, y_coord))
  R_alpha <- matrix(0, 2, 2)
  R_alpha[1, 1] <- cos(alpha)
  R_alpha[2, 2] <- cos(alpha)
  R_alpha[1, 2] <- -sin(alpha)
  R_alpha[2, 1] <- sin(alpha)
  new_coord <- R_alpha %*% coord
  return(list(x.coord = new_coord[1, ], y.coord = new_coord[2, ]))
}

create.ellispse.df = function(res, method){
  ncol <- 2
  nrow <- length(res$x.coord)
  df <- array(0, dim = c(nrow, 3))
  df[, 1] <- rep(method, nrow)
  df[, 2] <- res$x.coord
  df[, 3] <- res$y.coord
  return(as.data.frame(df))
}

#Outlier-free
true <- NULL
true$cov <- cov(xmat[-gt, ])
true$center <- apply(xmat[-gt, ], MARGIN = 2, FUN = function(h){mean(h, na.rm = TRUE)})

#OGK
obj.ogk <- CovOgk(xmat)
ogk.center <- getCenter(ogk)
ogk.cov <- getCov(ogk)

#my_ogk
my.ogk <- my_ogk(xmat)

df <- 

plot(xmat[, 1], xmat[, 2], pch = 19, cex = 0.5)
points(my_ogk$center[1], my_ogk$center[2], pch = 19, cex = 0.5, col = "red")
points(ogk.center[1], ogk.center[2], pch = 19, cex = 0.5, col = "green")
points(true.center[1], true.center[2], pch = 19, cex = 0.5, col = "blue")
points(new_coord_my_ogk[1, ], new_coord_my_ogk[2, ], col = "red", pch = 19, cex = 0.1)
points(xmat[gt, 1], xmat[gt, 2], pch = 19, cex = 0.5, col = "lightblue")
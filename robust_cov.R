mrob = function(h){tmp <- robustbase::scaleTau2(h, mu.too = TRUE); return(tmp[1])}
srob = function(h){robustbase::scaleTau2(h)}

ogk_R = function(x, m = mrob, s = srob){
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
      if (k != j){
        U[j, k] = 1/(4 * (s(ymat[, j] + ymat[, k])^2 - s(ymat[, j] - ymat[, k])^2))
      } else {
        U[j, k] = 1
      }
    }
  }
  E <- eigen(U, symmetric = TRUE)$vectors
  V <- ymat %*% E
  lambda <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p){
    lambda[i, i] <- s(V[, i])^2
  }
  mvec <- matrix(0, nrow = p, ncol = 1)
  for (i in 1:p){
    mvec[i, 1] <- mrob(V[, i])
  }
  mu_vec <- E %*% mvec
  sigmat <- E %*% lambda %*% t(E)
  
  return(list(center = D %*% mu_vec, cov = D %*% sigmat %*% t(D),U = U,V=V,E=E))
}

create.coord.df = function(x, method, ci = confidence, n.pts = 1000){
  if (method != "neutral" && method != "outlier"){
    xcenter <- x$center
    xcov <- x$cov
    p <- length(xcenter)
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
    df <- array(0, dim = c(n.pts, 3))
    df[, 1] <- rep(method, n.pts)
    df[, 2] <- new_coord[1, ]
    df[, 3] <- new_coord[2, ]
  } else {
    nrow <- nrow(x)
    df <- array(0, dim = c(nrow, 3))
    df[, 1] <- rep(method, nrow)
    df[, 2] <- x[, 1]
    df[, 3] <- x[, 2]
  }
  return(as.data.frame(df))
}

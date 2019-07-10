ishigami.fun <- function(X) {
  A <- 7
  B <- 0.1
  sin(X[, 1]) + A * sin(X[, 2]) ^ 2 + B * X[, 3] ^ 4 * sin(X[, 1])
}

X <- matrix(runif(3 * 20, -pi, pi), ncol = 3)
Y <- ishigami.fun(X)

Ydat <- Y
Xdat <- X

estimation <-
  topsa::topsa(
    Ydat = Y,
    Xdat = X,
    method = "VR",
    threshold.radius = 0.2,
    knearest = 25
  )

topsa::plot(estimation)
topsa::print(estimation)

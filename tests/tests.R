ishigami.fun <- function(X) {
A <- 7
B <- 0.1
sin(X[, 1]) + A * sin(X[, 2])^2 + B * X[, 3]^4 * sin(X[, 1])
}

a <- -pi
b <- pi
n <- 100
X1 <- runif(n, a, b)
X2 <- runif(n, a, b)
X3 <- runif(n, a, b)
Xdat <- data.frame(X1, X2, X3)
Ydat <- data.frame(Y = ishigami.fun(Xdat))


message("Ishigami Case")
ExIshigamiVR <- TopSA::TopSA(
  Ydat =  Ydat,
  Xdat = Xdat,
  dimension = 3,
  threshold = 0.01,
  knearest = 25,
  method = "VR",
  mc.cores = 1
)

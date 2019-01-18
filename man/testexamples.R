#Test for LatinR

a <- -1
b <- 1
n <- 500
X1 <- runif(n, a, b)
X2 <- runif(n, a, b)
Xr <- matrix(runif(n * 3, a, b), nrow = n)
X <- cbind(X1, X2, Xr)
Y <- 2 * X1 +  X2
# plot(X[, 1], Y)
# plot(X[, 2], Y)
# plot(X[, 3], Y)


radius <- rep(2, times = ncol(X))

message("Linear Case")
ExLinear <- TopSA::TopSA(
  Y = Y,
  X = X,
  radius = radius,
  dimension = 3,
  alpha = 0.05
)
g1 <- ExLinear$HOMOLOGY[[1]]
pdf(file = "man/exLinear-1.pdf")
plot.sobolManifold(g1, n=10000)
dev.off()
g2 <- ExLinear$HOMOLOGY[[2]]
pdf(file = "man/exLinear-2.pdf")
plot.sobolManifold(g2, n=10000)
dev.off()
g3 <- ExLinear$HOMOLOGY[[3]]
pdf(file = "man/exLinear-3.pdf")
plot.sobolManifold(g3, n=10000)
dev.off()


a <- -1
b <- 1
n <- 1000
X1 <- runif(n, a, b)
X2 <- runif(n, a, b)
X3 <- runif(n, a, b)
X <- cbind(X1, X2, X3)
Y <- X1 +  X2^4
# plot(X[,1],Y)
# plot(X[,2],Y)


radius <- rep(2, times = ncol(X))

message("Quartic Case")
ExQuartic <- sobolManifold(
  Y = Y,
  X = X,
  radius = radius,
  dimension = 3,
  alpha = 0.05
)


g1 <- ExQuartic$HOMOLOGY[[1]]
pdf(file = "man/exQuartic-1.pdf")
plot.sobolManifold(g1, n=10000)
dev.off()
g2 <- ExQuartic$HOMOLOGY[[2]]
pdf(file = "man/exQuartic-2.pdf")
plot.sobolManifold(g2, n=10000)
dev.off()
g3 <- ExQuartic$HOMOLOGY[[3]]
pdf(file = "man/exQuartic-3.pdf")
plot.sobolManifold(g3, n=10000)
dev.off()



a <- -1
b <- 1
n <- 500
theta <- runif(n, 0, 2*pi)
r <- (sqrt(runif(n))) * (0.5) + 0.5
X2 <- runif(n, a, b)
X1 <- r*cos(theta)
Y <- r*sin(theta)
X <- cbind(X1, X2)
# X <- scales::rescale(X, to = c(0, 1))
# Y <- scales::rescale(Y, to = c(0, 1))
# plot(X[,1],Y, asp = 1)
# plot(X[,2],Y)

radius <- c(0.5,1)

message("Circle 1 Case")
ExCirc1hole <- TopSA::TopSA(
  Y = Y,
  X = X,
  radius = radius,
  dimension = 3,
  alpha = 0.05
)

g1 <- ExCirc1hole$HOMOLOGY[[1]]
pdf(file = "man/exCircle-1.pdf")
plot.sobolManifold(g1, n=10000)
dev.off()

g2 <- ExCirc1hole$HOMOLOGY[[2]]
pdf(file = "man/exCircle-2.pdf")
plot.sobolManifold(g2, n=10000)
dev.off()

#
# alphaVec <- seq(0.05, 0.5, by = 0.05)
# ExCirc1holeAlpha <- NULL
#
# for (a in alphaVec) {
#   circ <- sobolManifold(
#     Y = Y,
#     X = X,
#     radius = radius,
#     dimension = 3,
#     alpha = a
#   )
#   ExCirc1holeAlpha <- cbind(a, circ$index[, 3])
# }





n <- 500
theta <- runif(n, 0, 2*pi)
r1 <- sqrt(runif(n)) + 1
r2 <- sqrt(runif(n)) * 0.5 + 0.5
X2 <- runif(2*n, -1, 1)

X1c1 <- r1*cos(theta)
Yc1 <- r1*sin(theta)
X1c2 <- r2*cos(theta)+1.5
Yc2 <- r2*sin(theta)+1.5

X1 <- c(X1c1, X1c2)
Y <- c(Yc1, Yc2)

X <- cbind(X1, X2)
# X <- scales::rescale(X, to = c(0, 1))
# Y <- scales::rescale(Y, to = c(0, 1))
# plot(X[,1],Y, asp=1)
# plot(X[,2],Y)

radius <- c(0.5,1)

message("Circle 2 Case")
ExCirc2holes <- TopSA::TopSA(
  Y = Y,
  X = X,
  radius = radius,
  dimension = 3,
  alpha = 0.05
)

g1 <- ExCirc2holes$HOMOLOGY[[1]]
pdf(file = "man/exCircle2holes-1.pdf")
plot.sobolManifold(g1, n=10000)
dev.off()

g2 <- ExCirc2holes$HOMOLOGY[[2]]
pdf(file = "man/exCircle2holes-2.pdf")
plot.sobolManifold(g2, n=10000)
dev.off()


a <- -pi
b <- pi
n <- 500
X1 <- runif(n, a, b)
X2 <- runif(n, a, b)
X3 <- runif(n, a, b)
X4 <- rnorm(n)
X <- cbind(X1, X2, X3, X4)
Y <- sensitivity::ishigami.fun(X)
# X <- scales::rescale(X, to = c(0,1))
# Y <- scales::rescale(Y, to = c(0,1))
# plot(scales  s::rescale(X[,1], to=c(0,1)),rescale(Y, to=c(0,1)), asp = 1)
# plot(scales::rescale(X[,2], to=c(0,1)),rescale(Y, to=c(0,1)), asp = 1)
# plot(scales::rescale(X[,3], to=c(0,1)),rescale(Y, to=c(0,1)), asp = 1)
# plot(X[,1],Y, asp = 0, xlim = c(-3.2, 3.2))
# plot(X[,2],Y, asp = 0, xlim = c(-3.2, 3.2))
# plot(X[,3],Y, asp = 0, xlim = c(-3.2, 3.2))

message("Ishigami Case")
ExIshigami <- TopSA::TopSA(
  Y = Y,
  X = X,
  radius = NULL,
  dimension = 3,
  alpha = 0.07
)
g1 <- ExIshigami$HOMOLOGY[[1]]
pdf(file = "man/exIshigami-1.pdf")
plot.sobolManifold(g1, n=10000)
dev.off()
g2 <- ExIshigami$HOMOLOGY[[2]]
pdf(file = "man/exIshigami-2.pdf")
plot.sobolManifold(g2, n=10000)
dev.off()
g3 <- ExIshigami$HOMOLOGY[[3]]
pdf(file = "man/exIshigami-3.pdf")
plot.sobolManifold(g3, n=10000)
dev.off()


save(file = "man/examples.Rdata",
     ExIshigami,
     ExCirc1hole,
     ExCirc2holes,
     ExLinear,
     ExQuartic)


# library(VineCopula)
# u <- pobs(as.matrix(cbind(X[,1],Y)))[,1]
# v <- pobs(as.matrix(cbind(X[,2],Y)))[,2]
# selectedCopula <- BiCopSelect(u,v,familyset=NA)
# selectedCopula



# X2[50] <- -1
# X <- cbind(X1,X2)
# Y <- X1 +  X2 ^ 2
# Y[50] <- 0.5
# plot(X1,Y)
# plot(X2,Y)
# plot(X3,Y)
#
# radius <- 0.2
# dimension <- 3
#
# aa <- sobolManifold(
#   Y = Y,
#   X = X,
#   radius = radius,
#   dimension = 3
# )
#
#
#
# g <- aa$HOMOLOGY[[1]]$graph
# plot.igraph(
#   mm,
#   rescale = FALSE,
#   axes = TRUE,
#   asp = 0,
#   ylim = range(V(g)$x),
#   xlim = range(V(g)$y),
#   vertex.size = 1,
#   vertex.label.cex = 1,
#   edge.arrow.size = 1,
#   edge.width = 0.5
# )
#
# g <- aa$HOMOLOGY[[2]]$graph
# plot.igraph(
#   g,
#   rescale = FALSE,
#   axes = TRUE,
#   asp = 0,
#   ylim = range(V(g)$x),
#   xlim = range(V(g)$y),
#   vertex.size = 1,
#   vertex.label.cex = 1,
#   edge.arrow.size = 1,
#   edge.width = 0.5
# )
#
# g <- aa$HOMOLOGY[[3]]$graph
# plot.igraph(
#   g,
#   rescale = FALSE,
#   axes = TRUE,
#   asp = 0,
#   ylim = range(V(g)$x),
#   xlim = range(V(g)$y),
#   vertex.size = 1,
#   vertex.label.cex = 1,
#   edge.arrow.size = 1,
#   edge.width = 0.5
# )


#Test for LatinR

a <- 0
b <- 5
n <- 100
X1 <- runif(n, a, b)
X2 <- runif(n, a, b)
X3 <- runif(n, a, b)
X4 <- runif(n, a, b)
X <- cbind(X1, X2, X3, X4)
Y <- X1 + 0.5 * X2

radius <- 0.5

Ex1Linear <- sobolManifold(
  Y = Y,
  X = X,
  radius = radius,
  dimension = 3
)






a <- -pi
b <- pi
n <- 1000
X1 <- runif(n, a, b)
X2 <- runif(n, a, b)
X3 <- runif(n, a, b)
X <- cbind(X1, X2, X3)
Y <- sensitivity::ishigami.fun(X)

Y <- sapply(1:ncol(X), function(i)
  return(Y[order(X[, i])]))    #order Y before X (or create a new variable for X)
X <- sapply(1:ncol(X), function(i)
  return(X[order(X[, i]), i]))


X2[50] <- -1
X <- cbind(X1,X2)
Y <- X1 +  X2 ^ 2
Y[50] <- 0.5
plot(X1,Y)
plot(X2,Y)
plot(X3,Y)

radius <- 0.2
dimension <- 3

aa <- sobolManifold(
  Y = Y,
  X = X,
  radius = radius,
  dimension = 3
)



g <- aa$HOMOLOGY[[1]]$graph
plot.igraph(
  mm,
  rescale = FALSE,
  axes = TRUE,
  asp = 0,
  ylim = range(V(g)$x),
  xlim = range(V(g)$y),
  vertex.size = 1,
  vertex.label.cex = 1,
  edge.arrow.size = 1,
  edge.width = 0.5
)

g <- aa$HOMOLOGY[[2]]$graph
plot.igraph(
  g,
  rescale = FALSE,
  axes = TRUE,
  asp = 0,
  ylim = range(V(g)$x),
  xlim = range(V(g)$y),
  vertex.size = 1,
  vertex.label.cex = 1,
  edge.arrow.size = 1,
  edge.width = 0.5
)

g <- aa$HOMOLOGY[[3]]$graph
plot.igraph(
  g,
  rescale = FALSE,
  axes = TRUE,
  asp = 0,
  ylim = range(V(g)$x),
  xlim = range(V(g)$y),
  vertex.size = 1,
  vertex.label.cex = 1,
  edge.arrow.size = 1,
  edge.width = 0.5
)

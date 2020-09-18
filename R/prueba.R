
ishigami.fun <- function(X) {
  A <- 7
  B <- 0.1
  sin(X[, 1]) + A * sin(X[, 2])^2 + B * X[, 3]^4 * sin(X[, 1])
}

X <- matrix(runif(3*1000, -pi, pi), ncol = 3)
Y <- ishigami.fun(X)

RotMat <- function(angle){
  matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), nrow=2, ncol=2)
}

library(TDA)
Xdat <- X
Ydat <- Y

X <- Xdat[,3]
Y <- Ydat
angle <- 0
X <- as.matrix(scales::rescale(X , to = c(0, 1)), ncol = 1)
Y <- as.matrix(scales::rescale(Y , to = c(0, 1)), ncol = 1)
meanX <- mean(X)
meanY <- mean(Y)
#Se centra para el calculo
XsYs <- data.frame(X-meanX,Y-meanY)
Rotated <- as.matrix(XsYs) %*% RotMat(-angle)
#Se descentra para el calculo
X <- Rotated[,1] + meanX
Y <- Rotated[,2] + meanY
a1 <- alphaComplexDiag(cbind(Y,X))


X <- Xdat[,3]
Y <- Ydat
angle <- 0
X <- as.matrix(scales::rescale(X , to = c(0, 1)), ncol = 1)
Y <- as.matrix(scales::rescale(Y , to = c(0, 1)), ncol = 1)
meanX <- mean(X)
meanY <- mean(Y)
#Se centra para el calculo
XsYs <- data.frame(X-meanX,Y-meanY)
Rotated <- as.matrix(XsYs) %*% RotMat(-angle)
#Se descentra para el calculo
X <- Rotated[,1] + meanX
Y <- Rotated[,2] + meanY
ix <- Y<1/4
Y <- Y[!ix]
X <- X[!ix]
a2 <- alphaComplexDiag(cbind(Y,X))

X <- Xdat[,3]
Y <- Ydat
angle <- 0
X <- as.matrix(scales::rescale(X , to = c(0, 1)), ncol = 1)
Y <- as.matrix(scales::rescale(Y , to = c(0, 1)), ncol = 1)
meanX <- mean(X)
meanY <- mean(Y)
#Se centra para el calculo
XsYs <- data.frame(X-meanX,Y-meanY)
Rotated <- as.matrix(XsYs) %*% RotMat(-angle)
#Se descentra para el calculo
X <- Rotated[,1] + meanX
Y <- Rotated[,2] + meanY
ix <- Y<1/2
Y <- Y[!ix]
X <- X[!ix]
a3 <- alphaComplexDiag(cbind(Y,X))

X <- Xdat[,3]
Y <- Ydat
angle <- 0
X <- as.matrix(scales::rescale(X , to = c(0, 1)), ncol = 1)
Y <- as.matrix(scales::rescale(Y , to = c(0, 1)), ncol = 1)
meanX <- mean(X)
meanY <- mean(Y)
#Se centra para el calculo
XsYs <- data.frame(X-meanX,Y-meanY)
Rotated <- as.matrix(XsYs) %*% RotMat(-angle)
#Se descentra para el calculo
X <- Rotated[,1] + meanX
Y <- Rotated[,2] + meanY
ix <- Y<3/4
Y <- Y[!ix]
X <- X[!ix]
a4 <- alphaComplexDiag(cbind(Y,X))


plot(X,Y)
plot.diagram(a1[["diagram"]])
plot.diagram(a2[["diagram"]])
plot.diagram(a3[["diagram"]])
plot.diagram(a4[["diagram"]])

bottleneck(a1[["diagram"]],a2[["diagram"]])
bottleneck(a2[["diagram"]],a3[["diagram"]])
bottleneck(a3[["diagram"]],a4[["diagram"]])


particion <- seq(0,1, length.out = 20)
a<- list()

for(k in particion){
X <- Xdat[,2]
Y <- Ydat
angle <- 0
X <- as.matrix(scales::rescale(X , to = c(0, 1)), ncol = 1)
Y <- as.matrix(scales::rescale(Y , to = c(0, 1)), ncol = 1)
meanX <- mean(X)
meanY <- mean(Y)
#Se centra para el calculo
XsYs <- data.frame(X-meanX,Y-meanY)
Rotated <- as.matrix(XsYs) %*% RotMat(-angle)
#Se descentra para el calculo
X <- Rotated[,1] + meanX
Y <- Rotated[,2] + meanY
ix <- Y<k
Y <- Y[!ix]
X <- X[!ix]
a <- c(a,alphaComplexDiag(cbind(Y,X)))
}

b<- numeric()
for(k in 1:(length(a)-1)){

 b[k] <-bottleneck(a[k][["diagram"]],a[k+1][["diagram"]])
}


k=particion[15]

plot(b)

plot(X,Y)



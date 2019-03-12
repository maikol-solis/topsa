library(TopSA)

#Test for LatinR

a <- -1
b <- 1
n <- 1000
X1 <- runif(n, a, b)
X2 <- runif(n, a, b)
Xnoise <- matrix(runif(n * 3, a, b), nrow = n)
Xdat <- data.frame(X1, X2, Xnoise = Xnoise)
Ydat <-data.frame(Y = 2 * X1 +  X2)
# plot(X[, 1], Y)
# plot(X[, 2], Y)
# plot(X[, 3], Y)


message("Linear Case")
ExLinearVR <- TopSA::TopSA(
  Ydat = Ydat,
  Xdat = Xdat,
  dimension = 3,
  threshold = 0.05,
  method = "VR"
)

pdf(file = "man/exLinearVR-reflection-1.pdf")
plot(ExLinearVR,nvar = 1,with.reflection = TRUE, asp=0)
dev.off()
pdf(file = "man/exLinearVR-sym-diff-1.pdf")
plot(ExLinearVR,nvar = 1,symmetric.diff = TRUE, asp=0)
dev.off()
pdf(file = "man/exLinearVR-reflection-2.pdf")
plot(ExLinearVR,nvar = 2, with.reflection = TRUE, asp=0)
dev.off()
pdf(file = "man/exLinearVR-sym-diff-2.pdf")
plot(ExLinearVR,nvar = 2, symmetric.diff = TRUE, asp=0)
dev.off()
pdf(file = "man/exLinearVR-reflection-3.pdf")
plot(ExLinearVR,nvar = 3, with.reflection = TRUE, asp=0)
dev.off()
pdf(file = "man/exLinearVR-sym-diff-3.pdf")
plot(ExLinearVR,nvar = 3, symmetric.diff = TRUE, asp=0)
dev.off()

LinearVR <- print(ExLinearVR, only.return.table = T)
xtable::print.xtable(xtable::xtable(LinearVR, align = c("c","c","c","c","c","c"),file = "man/exLinearVR-table.tex"))


message("Linear Case")
ExLineardelanauy <- TopSA::TopSA(
  Ydat = Ydat,
  Xdat = Xdat,
  dimension = 3,
  threshold = 0.01,
  method = "delanauy"
)

pdf(file = "man/exLineardelanauy-reflection-1.pdf")
plot(ExLineardelanauy,nvar = 1,with.reflection = TRUE, asp=0)
dev.off()
pdf(file = "man/exLineardelanauy-sym-diff-1.pdf")
plot(ExLineardelanauy,nvar = 1,symmetric.diff = TRUE, asp=0)
dev.off()
pdf(file = "man/exLineardelanauy-reflection-2.pdf")
plot(ExLineardelanauy,nvar = 2, with.reflection = TRUE, asp=0)
dev.off()
pdf(file = "man/exLineardelanauy-sym-diff-2.pdf")
plot(ExLineardelanauy,nvar = 2, symmetric.diff = TRUE, asp=0)
dev.off()
pdf(file = "man/exLineardelanauy-reflection-3.pdf")
plot(ExLineardelanauy,nvar = 3, with.reflection = TRUE, asp=0)
dev.off()
pdf(file = "man/exLineardelanauy-sym-diff-3.pdf")
plot(ExLineardelanauy,nvar = 3, symmetric.diff = TRUE, asp=0)
dev.off()

Lineardelanauy <- print(ExLineardelanauy, only.return.table = T)
xtable::print.xtable(xtable::xtable(Lineardelanauy, align = c("c","c","c","c","c","c"),file = "man/exLineardelanauy-table.tex"))


# a <- -1
# b <- 1
# n <- 1000
# X1 <- runif(n, a, b)
# X2 <- runif(n, a, b)
# X3 <- runif(n, a, b)
# X <- cbind(X1, X2, X3)
# Y <- X1 +  X2^4
# # plot(X[,1],Y)
# # plot(X[,2],Y)
#
#
# radius <- rep(2, times = ncol(X))
#
# message("Quartic Case")
# ExQuartic <- sobolManifold(
#   Y = Y,
#   X = X,
#   radius = radius,
#   dimension = 3,
#   alpha = 0.05
# )
#
#
# g1 <- ExQuartic$HOMOLOGY[[1]]
# pdf(file = "man/exQuartic-1.pdf")
# plot.sobolManifold(g1, n=10000)
# dev.off()
# g2 <- ExQuartic$HOMOLOGY[[2]]
# pdf(file = "man/exQuartic-2.pdf")
# plot.sobolManifold(g2, n=10000)
# dev.off()
# g3 <- ExQuartic$HOMOLOGY[[3]]
# pdf(file = "man/exQuartic-3.pdf")
# plot.sobolManifold(g3, n=10000)
# dev.off()




a <- -1
b <- 1
n <- 1000
theta <- runif(n, 0, 2*pi)
r <- (sqrt(runif(n))) * (0.5) + 0.5
X2 <- runif(n, a, b)
X1 <- r*cos(theta)
Ydat <- data.frame(Y = r*sin(theta))
Xdat <- data.frame(X1, X2)


message("Circle 1 Case")
ExCirc1holeVR <- TopSA::TopSA(
  Ydat = Ydat,
  Xdat = Xdat,
  dimension = 3,
  threshold = 0.05,
  knearest = 20,
  method = "VR"
)

pdf(file = "man/ExCirc1holeVR-reflection-1.pdf")
plot(ExCirc1holeVR,nvar = 1,with.reflection = TRUE, asp = 1)
dev.off()
pdf(file = "man/ExCirc1holeVR-sym-diff-1.pdf")
plot(ExCirc1holeVR,nvar = 1,symmetric.diff = TRUE)
dev.off()
pdf(file = "man/ExCirc1holeVR-reflection-2.pdf")
plot(ExCirc1holeVR,nvar = 2, with.reflection = TRUE)
dev.off()
pdf(file = "man/ExCirc1holeVR-sym-diff-2.pdf")
plot(ExCirc1holeVR,nvar = 2, symmetric.diff = TRUE)
dev.off()

Circle1holeVR <- print(ExCirc1holeVR, only.return.table = T)
xtable::print.xtable(xtable::xtable(Circle1holeVR,align = c("c","c","c","c","c","c"),file = "man/exCircle1holeVR-table.tex"))

ExCirc1holedelanauy <- TopSA::TopSA(
  Ydat = Ydat,
  Xdat = Xdat,
  dimension = 3,
  threshold = 0.01,
  method = "delanauy"
)

pdf(file = "man/ExCirc1holedelanauy-reflection-1.pdf")
plot(ExCirc1holedelanauy,nvar = 1,with.reflection = TRUE, asp = 1)
dev.off()
pdf(file = "man/ExCirc1holedelanauy-sym-diff-1.pdf")
plot(ExCirc1holedelanauy,nvar = 1,symmetric.diff = TRUE)
dev.off()
pdf(file = "man/ExCirc1holedelanauy-reflection-2.pdf")
plot(ExCirc1holedelanauy,nvar = 2, with.reflection = TRUE)
dev.off()
pdf(file = "man/ExCirc1holedelanauy-sym-diff-2.pdf")
plot(ExCirc1holedelanauy,nvar = 2, symmetric.diff = TRUE)
dev.off()

Circle1holedelaunay <- print(ExCirc1holedelanauy, only.return.table = T)
xtable::print.xtable(xtable::xtable(Circle1holedelaunay,align = c("c","c","c","c","c","c"),file = "man/exCircle1holedelaunay-table.tex"))



n <- 1000
theta <- runif(n, 0, 2*pi)
r1 <- sqrt(runif(n)) + 1
r2 <- sqrt(runif(n)) * 0.5 + 0.5
X2 <- runif(2*n, -1, 1)

X1c1 <- r1*cos(theta)
Yc1 <- r1*sin(theta)
X1c2 <- r2*cos(theta)+1.5
Yc2 <- r2*sin(theta)+1.5

X1 <- c(X1c1, X1c2)
Ydat <- data.frame(Y = c(Yc1, Yc2))

Xdat <- data.frame(X1, X2)

message("Circle 2 Case")
ExCirc2holesVR <- TopSA::TopSA(
  Ydat = Ydat,
  Xdat = Xdat,
  dimension = 3,
  threshold = 0.05,
  knearest = 20,
  method = "VR"
)

pdf(file = "man/ExCirc2holesVR-reflection-1.pdf")
plot(ExCirc2holesVR,nvar = 1,with.reflection = TRUE, asp = 1)
dev.off()
pdf(file = "man/ExCirc2holesVR-sym-diff-1.pdf")
plot(ExCirc2holesVR,nvar = 1,symmetric.diff = TRUE)
dev.off()
pdf(file = "man/ExCirc2holesVR-reflection-2.pdf")
plot(ExCirc2holesVR,nvar = 2, with.reflection = TRUE)
dev.off()
pdf(file = "man/ExCirc2holesVR-sym-diff-2.pdf")
plot(ExCirc2holesVR,nvar = 2, symmetric.diff = TRUE)
dev.off()



ExCirc2holesdelanauy <- TopSA::TopSA(
  Ydat = Ydat,
  Xdat = Xdat,
  dimension = 3,
  threshold = 0.005,
  method = "delanauy"
)

pdf(file = "man/ExCirc2holesdelanauy-reflection-1.pdf")
plot(ExCirc2holesdelanauy,nvar = 1,with.reflection = TRUE, asp = 1)
dev.off()
pdf(file = "man/ExCirc2holesdelanauy-sym-diff-1.pdf")
plot(ExCirc2holesdelanauy,nvar = 1,symmetric.diff = TRUE)
dev.off()
pdf(file = "man/ExCirc2holesdelanauy-reflection-2.pdf")
plot(ExCirc2holesdelanauy,nvar = 2, with.reflection = TRUE)
dev.off()
pdf(file = "man/ExCirc2holesdelanauy-sym-diff-2.pdf")
plot(ExCirc2holesdelanauy,nvar = 2, symmetric.diff = TRUE)
dev.off()



a <- -pi
b <- pi
n <- 1000
X1 <- runif(n, a, b)
X2 <- runif(n, a, b)
X3 <- runif(n, a, b)
Xdat <- data.frame(X1, X2, X3)
Ydat <- data.frame(Y = sensitivity::ishigami.fun(Xdat))


message("Ishigami Case")
ExIshigamiVR <- TopSA::TopSA(
  Ydat =  Ydat,
  Xdat = Xdat,
  dimension = 3,
  threshold = 0.03,
  knearest = 25,
  method = "VR"
)

pdf(file = "man/ExIshigamiVR-reflection-1.pdf")
plot(ExIshigamiVR, nvar = 1, with.reflection = TRUE, asp=0)
dev.off()
pdf(file = "man/ExIshigamiVR-sym-diff-1.pdf")
plot(ExIshigamiVR, nvar = 1, symmetric.diff = TRUE, asp=0)
dev.off()
pdf(file = "man/ExIshigamiVR-reflection-2.pdf")
plot(ExIshigamiVR, nvar = 2, with.reflection = TRUE, asp=0)
dev.off()
pdf(file = "man/ExIshigamiVR-sym-diff-2.pdf")
plot(ExIshigamiVR, nvar = 2, symmetric.diff = TRUE, asp=0)
dev.off()
pdf(file = "man/ExIshigamiVR-reflection-3.pdf")
plot(ExIshigamiVR, nvar = 3, with.reflection = TRUE, asp=0)
dev.off()
pdf(file = "man/ExIshigamiVR-sym-diff-3.pdf")
plot(ExIshigamiVR, nvar = 3, symmetric.diff = TRUE, asp=0)
dev.off()
pdf(file = "man/ExIshigamiVR-reflection-4.pdf")
plot(ExIshigamiVR, nvar = 4, with.reflection = TRUE, asp=0)
dev.off()
pdf(file = "man/ExIshigamiVR-sym-diff-4.pdf")
plot(ExIshigamiVR, nvar = 4, symmetric.diff = TRUE, asp=0)
dev.off()

IshigamiVR <- print(ExIshigamiVR, only.return.table = T)
xtable::print.xtable(xtable::xtable(IshigamiVR,align = c("c","c","c","c","c","c"),file = "man/exIshigamiVR-table.tex"))



ExIshigamidelanauy <- TopSA::TopSA(
  Ydat =  Ydat,
  Xdat = Xdat,
  dimension = 3,
  threshold = 0.02,
  method = "delanauy"
)

pdf(file = "man/ExIshigamidelanauy-reflection-1.pdf")
plot(ExIshigamidelanauy, nvar = 1, with.reflection = TRUE, asp=0)
dev.off()
pdf(file = "man/ExIshigamidelanauy-sym-diff-1.pdf")
plot(ExIshigamidelanauy, nvar = 1, symmetric.diff = TRUE, asp=0)
dev.off()
pdf(file = "man/ExIshigamidelanauy-reflection-2.pdf")
plot(ExIshigamidelanauy, nvar = 2, with.reflection = TRUE, asp=0)
dev.off()
pdf(file = "man/ExIshigamidelanauy-sym-diff-2.pdf")
plot(ExIshigamidelanauy, nvar = 2, symmetric.diff = TRUE, asp=0)
dev.off()
pdf(file = "man/ExIshigamidelanauy-reflection-3.pdf")
plot(ExIshigamidelanauy, nvar = 3, with.reflection = TRUE, asp=0)
dev.off()
pdf(file = "man/ExIshigamidelanauy-sym-diff-3.pdf")
plot(ExIshigamidelanauy, nvar = 3, symmetric.diff = TRUE, asp=0)
dev.off()
pdf(file = "man/ExIshigamidelanauy-reflection-4.pdf")
plot(ExIshigamidelanauy, nvar = 4, with.reflection = TRUE, asp=0)
dev.off()
pdf(file = "man/ExIshigamidelanauy-sym-diff-4.pdf")
plot(ExIshigamidelanauy, nvar = 4, symmetric.diff = TRUE, asp=0)
dev.off()

Ishigamidelaunay <- print(ExIshigamidelanauy, only.return.table = T)
xtable::print.xtable(xtable::xtable(Ishigamidelaunay,align = c("c","c","c","c","c","c"),file = "man/exIshigamidelaunay-table.tex"))


save(
  file = "man/examples.Rdata",
  ExCirc1holedelanauy,
  ExCirc1holeVR,
  ExCirc2holesdelanauy,
  ExCirc2holesVR,
  ExIshigamidelanauy,
  ExIshigamiVR,
  ExLineardelanauy,
  ExLinearVR
)



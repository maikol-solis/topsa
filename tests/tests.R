# library(topsa)
# library(ggplot2)
#
# # Linear example
#
# linear.fun <- function(X) {
#   - 2 * X[, 1] +  X[, 2]
# }
#
# X <- matrix(runif(3 * 1000, -1, 1), ncol = 3)
# Y <- linear.fun(X)
#
# estimation <-
#   topsa::topsa(
#     Ydat = Y,
#     Xdat = X,
#     method = "VR",
#     threshold.radius = 0.05,
#     knearest = 10
#   )
#
# p <- topsa::plot(estimation)
# p[[1]]
# p[[2]]
# p[[3]]
# topsa::print(estimation)
# un

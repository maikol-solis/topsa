test_that("plot_topsa works", {
  linear.fun <- function(X) {
  -2 * X[, 1] +  X[, 2]
}

X <- data.frame(
  X1 = runif(20, -1, 1),
  X2 = runif(20, -1, 1),
  X3 = runif(20, -1, 1)
)
Y <- data.frame(linear.fun(X))

topsaObj <- topsa(Ydat = Y, Xdat = X)

p <- plot_topsa(topsaObj)

expect_s3_class(p, "ggplot")
})



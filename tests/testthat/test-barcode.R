test_that("barcode_plotter works", {
  linear.fun <- function(X) {
  -2 * X[, 1] +  X[, 2]
}

X <- data.frame(
  X1 = runif(20, -1, 1),
  X2 = runif(20, -1, 1),
  X3 = runif(20, -1, 1)
)
Y <- data.frame(linear.fun(X))

p <- barcode_plotter(Ydat = Y,
                    Xdat = X,
                    maxscale = 0.1,
                    mc.cores = 2)

expect_s3_class(p, "ggplot")
})

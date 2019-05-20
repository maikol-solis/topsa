#' Plot method for objects \code{TopSA}
#'
#' Plot the Complex, symmetric reflection and symmetric diffrence for and object
#' of class \code{\ĺink{TopSA}}
#'
#' @param TopSAObj an object of class \code{TopSA}
#' @param only.return.table  if \code{TRUE}, returns a data frame with the
#'   estimated values. Otherwise, print the data frame in console. Defaults to
#'   \code{FALSE}
#' @param ... further arguments passed to the \code{plot} function
#'
#' @return A plot of generated with the output of \code{\link{TopSA}} . For each
#' variable in the model, it creates the plot of the correponding complex, its
#' symmetric reflection and its symmetric difference.
#' @export
#'
#' @examples
#'
#' ishigami.fun <- function(X) {
#' A <- 7
#' B <- 0.1
#' sin(X[, 1]) + A * sin(X[, 2])^2 + B * X[, 3]^4 * sin(X[, 1])
#' }
#'
#' X <- matrix(runif(3*100, -pi, pi), ncol = 3)
#' Y <- ishigami.fun(X)
#'
#' estimation <- sobolnp(Y = Y, X = X, nboot = 5)
#'
#' plot(estimation)
#'
# @importFrom
#'

print <- function(TopSAObj, ...) {
  UseMethod("print", TopSAObj)
}

#' @export
#' @rdname print
print.TopSA <- function(TopSAObj, only.return.table = FALSE, ...) {
  sensitivity_table <- t(sapply(TopSAObj$results, function(x) {
    unlist(x[c(
      "threshold",
      "Manifold.Area",
      "Box.Area",
      "Geometric.Correlation",
      "Symmetric.Index"
    )])
  }))

  colnames(sensitivity_table) <-
    c('Threshold',
      'Manifold Area',
      'Box Area' ,
      'Geometric correlation',
      'Symmetric index')
  rownames(sensitivity_table) <- colnames(TopSAObj$Xdat)

  if(only.return.table == TRUE){
    return(sensitivity_table)
  }


  cat("\nCall:\n", deparse(TopSAObj[['call']]), "\n", sep = "")
  cat("\nMethod used:",deparse(TopSAObj[["call"]]$method), sep = "")
  cat("\nNumber of variables:", ncol(TopSAObj[['Xdat']]), "\n")
  cat("\nNumber of observations:", nrow(TopSAObj[['Ydat']]), "\n")
  cat("\nFirst order indices\n")
  print(sensitivity_table[, ])
  # return(sensitivity_table)
}

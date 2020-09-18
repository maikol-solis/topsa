#' print \code{topsa} objects
#'
#' Print method for objects of class \code{topsa}.
#'
#' @param topsaObj an object of class \code{topsa}
#' @param only.return.table  if \code{TRUE}, returns a data frame with the
#'   estimated values. Otherwise, print the data frame in console. Defaults to
#'   \code{FALSE}
#' @param ... further arguments passed to the \code{plot} function
#'
#' @return Print the threshold used, the box area, manifold embedding area, geometric
#' correlation index and symmetric sensitivity index for and object of class
#' \code{topsa}.
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
#' X <- matrix(runif(3*50, -pi, pi), ncol = 3)
#' Y <- ishigami.fun(X)
#'
#' estimation <- topsa(Ydat = Y, Xdat = X)
#'
#' print(estimation)

print_topsa <- function(topsaObj, only.return.table = FALSE, ...) {
  sensitivity_table <- t(sapply(topsaObj$results, function(x) {
    unlist(x[c(
      "threshold",
      "Manifold.Area",
      "Box.Area",
      "Geometric.R2",
      "Geometric.Sensitivity"
    )])
  }))

  colnames(sensitivity_table) <-
    c('Threshold',
      'Manifold Area',
      'Box Area' ,
      'Geometric R2',
      'Geometric Sensitivity')
  rownames(sensitivity_table) <- colnames(topsaObj$Xdat)

  if (only.return.table == TRUE) {
    return(sensitivity_table)
  }


  cat("\nCall:\n", deparse(topsaObj[['call']]), "\n", sep = "")
  cat("\nMethod used:", deparse(topsaObj[["call"]]$method), sep = "")
  cat("\nNumber of variables:", ncol(topsaObj[['Xdat']]), "\n")
  cat("\nNumber of observations:", nrow(topsaObj[['Ydat']]), "\n")
  cat("\nFirst order indices\n")
  print(sensitivity_table[,])
  # return(sensitivity_table)
}

#' Topological Sensitivity Analysis
#' @param Ydat numeric vector of responses in a model.
#' @param Xdat numeric matrix or data.frame of covariables.
# @param dimension number of homology groups to estimate. The only value accepted is \code{dimension = 3}, but in future release this will change.
#  @param threshold.area (no used for the moment) percent of the triangles bigger areas to keep.
#' @param threshold.radius percent of radius or sizes of triangles to keep.
#'   the homology complex. Defaults to `0.05`.
# @param knearest number of the nearest neigbors keep to built the homology
#'   complex
#' @param method type of method to build the homology complex. Two choices are
#'   accepted: \code{Alpha} o \code{VR} (Vietoris-Rips).
#'
#' @return A list of class \code{topsa} with the following elements:
#'
#' \describe{
#' \item{\strong{call}}{The function call.}
#' \item{\strong{Xdat}}{\code{X} input.}
#' \item{\strong{Ydat}}{\code{Y} output.}
#' \item{\strong{dimension}}{dimension to estimate the homology order.}
#' \item{\strong{threshold}}{cutoff level for the radius or area.}
#' \item{\strong{results}}{A list for each variable with:
#' \describe{
#' \item{\strong{threshold}}{threshold used to limit the radius or area.}
#' \item{\strong{Manifold_Area}}{geometrical area of the estimated manifold.}
#' \item{\strong{Box.Area}}{geometrical area of the estimated manifold.}
#' \item{\strong{Geometric.R2}}{geometric correlation between each
#' \code{x} and \code{y}.}
#' \item{\strong{Geometric.Sensitivity}}{symmetric sensitivity
#' index of each estimated manifold.}
#' \item{\strong{manifold_plot}}{a \code{sf}
#' object with the estimated manifold.}
#' } } }
#' @export
#'
#' @examples
#'
#'ishigami.fun <- function(X) {
#' A <- 7
#' B <- 0.1
#' sin(X[, 1]) + A * sin(X[, 2])^2 + B * X[, 3]^4 * sin(X[, 1])
#' }

#' X <- matrix(runif(3*50, -pi, pi), ncol = 3)
#' Y <- ishigami.fun(X)
#'estimation <- topsa(Ydat = Y, Xdat = X,method = "Alpha")
#' @importFrom methods is


topsa <-
  function(Ydat,
           Xdat,
           threshold.radius = rep(0.02, ncol(Xdat)),
           method = "Alpha",
           mc.cores = 2,
           angle = 0) {
    Xdat <- as.data.frame(Xdat)
    Ydat <- as.data.frame(Ydat)

    ANS <- list()
    ANS[['call']] <- match.call()
    ANS[['Xdat']] <- Xdat
    ANS[['Ydat']] <- Ydat

    # Xr <- matrix()
    # Yr <- matrix()
    # l <- lapply(seq_along(Xdat), function(k) {
    #   scales::rescale(cbind(Xdat[, k], Ydat[, 1]))
    # })
    #
    # lx <- lapply(l, function(x)
    #   x[, 1])
    #
    # Xr <- as.data.frame(do.call("cbind", lx))
    # Yr <- as.data.frame(sapply(Ydat, scales::rescale))
    # ANS[['Xr']] <- Xr
    # ANS[['Yr']] <- Yr
    ANS[['Xr']] <- as.data.frame(lapply(Xdat, scales::rescale))
    ANS[['Yr']] <- as.data.frame(lapply(Ydat, scales::rescale))
    ANS[['angle']] <- angle

    # ANS[['dimension']] <- dimension

    # message("Estimating persistance homology for each variable")
    # pb <-
    #   utils::txtProgressBar(min = 0,
    #                         max = ncol(Xdat),
    #                         style = 3)
    #
    # ANS[['persistence homology']] <-
    #   lapply(
    #     X = 1:ncol(Xdat),
    #     FUN =  function(i) {
    #       h <-
    #         TDAstats::calculate_homology(apply(cbind(Xdat[, i], Ydat), 2, scales::rescale), dim = 1)
    #       utils::setTxtProgressBar(pb, i)
    #       return(h)
    #     }
    #   )
    # close(pb)


    # lapply(
    #   X = 1:ncol(Xdat),
    #   FUN = function(i, Xdat, Ydat) {
    #     h <-
    #       TDAstats::calculate_homology(apply(cbind(Xdat[, i], Ydat), 2, scales::rescale))
    #   },
    #   Xdat = Xdat,
    #   Ydat = Ydat
    # )

    # ANS[['threshold']] <-
    #   if (method == "VR") {
    #     if (threshold.radius == -1) {
    #       threshold.radius <-  sapply(
    #         X = ANS[['persistence homology']],
    #         FUN = function(h) {
    #           # h <-  h[h[, "dimension"] == 1,]
    #
    #           diff_birth_death <-
    #             apply(h[, c("birth", "death")], 1, function(x) {
    #               x[2] - x[1]
    #             })
    #
    #           # q_diff <- quantile(diff_birth_death, 0.99)
    #           #
    #           # h <- as.matrix(h[diff_birth_death >= q_diff,])
    #           #
    #           # q_death <- quantile(h[, "death"], 0.99)
    #           #
    #           # h <- as.matrix(h[h[, "death"] >= q_death,])
    #
    #
    #
    #           order_diff <-
    #             order(diff_birth_death, decreasing = TRUE)
    #
    #           order_death <- order(h[,"death"],-h[, "death"], decreasing = TRUE)
    #
    #           order_in_both_seq <-
    #             cumsum(order_death == order_diff)
    #
    #           order_in_perfect_case <- 1:length(order_death)
    #
    #           idx_features <-
    #             order_death[order_in_both_seq == order_in_perfect_case]
    #
    #
    #           # Test if there are prominent features If lenght(idx_features)==0
    #           # then almost all the persistance is noise. In this case just take
    #           # the largest radius
    #           # if (length(idx_features) == 0) {
    #           #   threshold.radius <- max(h[, "death"])
    #           #
    #           # } else{
    #           # If not, then try to search the intersection of all the select
    #           # features. Once found this intersection, take the middle point.
    #           # intervals <- h[idx_features, c("birth", "death")]
    #
    #           # If intervals is a single vector, only assign l and r to the
    #           # first and second elements
    #           if (is.null(dim(h))) {
    #             intervals <- h[c("birth", "death")]
    #             l <- intervals[1]
    #             r <- intervals[2]
    #           } else{
    #             # Otherwise, try to find their intersection
    #             intervals <- h[, c("birth", "death")]
    #             l <- intervals[1, 1]
    #             r <- intervals[1, 2]
    #             # Find the intersection over all the intervals
    #             for (i in seq_len(nrow(intervals))) {
    #               # If no intersection exists
    #               if (intervals[i, 1] > r | intervals[i, 2] < l) {
    #                 l <- min(intervals[, 1])
    #                 r <- max(intervals[, 2])
    #                 break()
    #               } else{
    #                 # Else update the intersection
    #                 l <- max(l, intervals[i, 1])
    #                 r <- min(r, intervals[i, 2])
    #               }
    #
    #             }
    #
    #           }
    #           threshold.radius <- (l + r) / 2
    #
    #
    #           # }#end-else-if (length(idx_features) == 0)
    #
    #           return(threshold.radius)
    #
    #         }#end-FUN = function(i, Xdat, Ydat)
    #       )
    #     }#end-if-threshold.radius=-1
    #
    #         if (length(threshold.radius) == 1) {
    #           threshold.radius <- rep(threshold.radius, ncol(Xdat))
    #         } else if (length(threshold.radius) == ncol(Xdat)) {
    #           threshold.radius
    #         } else{
    #           stop("Please provide a numeric threshold vector of size 1 or ncol(Xdat)")
    #         }
    #
    #       }#end-if-method==VR


    if (length(threshold.radius) == 1) {
      threshold.radius <- rep(threshold.radius, ncol(Xdat))
    } else if (length(threshold.radius) < ncol(Xdat)) {
      stop("Please provide a numeric threshold vector of size 1 or ncol(Xdat)")
    }


    par.names = colnames(Xdat)	#gets parameters names
    if (is.null(colnames(Xdat))) {
      par.names = paste0('X', 1:ncol(Xdat))
    }


    message("Index estimation")


    # sensitivity_results <- try(parallel::mclapply(
    #   X = 1:ncol(Xdat),
    #   FUN = estimate_sensitivity_index,
    #   Ydat = Ydat,
    #   Xdat = Xdat,
    #   threshold = threshold.radius,
    #   method = method,
    #   mc.cores = mc.cores,
    #   angle = angle),
    # silent = T)



    # if (is(sensitivity_results, 'try-error')) {
      sensitivity_results <-
        lapply(
          X = 1:ncol(Xdat),
          FUN = estimate_sensitivity_index,
          Ydat = Ydat,
          Xdat = Xdat,
          threshold = threshold.radius,
          method = method,
          angle = angle
        )
    # }


    for (i in 1:ncol(Xdat)) {
      sensitivity_results[[i]]$xname <- colnames(Xdat)[i]
      sensitivity_results[[i]]$yname <- colnames(Ydat)
    }

    ANS[["results"]] <- sensitivity_results
    class(ANS) <- 'topsa'
    return(ANS)


  }

estimate_sensitivity_index <- function(ivar,
                                       Ydat,
                                       Xdat,
                                       dimension,
                                       # knearest,
                                       threshold,
                                       method,
                                       angle) {
  constructHOMOLOGY <-
    function (ivar,
              Ydat,
              Xdat,
              dimension,
              threshold,
              method,
              angle) {
      Y <- as.matrix(Ydat)
      X <- as.matrix(Xdat[,ivar])#CAMBIO

      idx <- order(X, Y)
      X <- X[idx, ]
      Y <- Y[idx, ]

      YX <- scales::rescale(cbind(Y, X))
      Y <- YX[, 1]
      X <- YX[, 2]
      # X <- as.matrix(scales::rescale(X , to = c(0, 1)), ncol = 1)
      # Y <- as.matrix(scales::rescale(Y , to = c(0, 1)), ncol = 1)

      # meanX <- mean(X)
      # meanY <- mean(Y)
      # #Se centra para el calculo
      # XsYs <- data.frame(X - meanX, Y - meanY)
      # Rotated <- as.matrix(XsYs) %*% RotMat(-angle)
      # #Se descentra para el calculo
      # X <- Rotated[, 1] + meanX
      # Y <- Rotated[, 2] + meanY
      #
      # YX <- scales::rescale(cbind(Y, X))
      # Y <- YX[, 1]
      # X <- YX[, 2]



      if (method == "Alpha") {
        Filtration <-
          TDA::alphaComplexFiltration(cbind(Y, X), printProgress = FALSE)
        cmplx <-
          Filtration$cmplx[Filtration$values <= threshold[ivar]]
      } else if (method == "VR") {
        Filtration <- TDA::ripsFiltration(
          cbind(Y, X),
          maxdimension = 1,
          maxscale = threshold[ivar],
          printProgress = FALSE
        )
        cmplx <- Filtration$cmplx
      } else{
        Filtration <- NULL
        stop("No method defined")
      }


      idx_triangles <- lengths(Filtration$cmplx) == 3
      clq <- cmplx[idx_triangles]
      # clq <- igraph::cliques(graphBase, min = 3, max = 3)
      clq <- matrix(unlist(clq), ncol = 3, byrow = TRUE)
      clq <- cbind(clq, clq[, 1])
      clq <- clq[order(clq[, 1], -clq[, 2], clq[, 3]),]

      clq_polygons <-
        lapply(
          X = seq_len(nrow(clq)),
          FUN =  function(i) {
            p <- data.frame(id = i, x = X[clq[i, ]], y = Y[clq[i, ]])
          }
        )

      clq_polygons <- do.call("rbind", clq_polygons)

      clq_polygons <-
        sfheaders::sfc_polygon(clq_polygons, polygon_id = "id")



      clq_polygons <- sf::st_union(clq_polygons)
      clq_polygons <- sf::st_cast(clq_polygons, "POLYGON")
      clq_polygons <- clq_polygons[sf::st_area(clq_polygons) > 0.01]
      clq_polygons <- sf::st_buffer(clq_polygons, -0.01)
      clq_polygons <- sf::st_buffer(clq_polygons, 0.01)
      clq_polygons <- sf::st_union(clq_polygons)



      return(
        list(
          manifold_unioned = clq_polygons,
          # neigborhood.distance = neigborhood.distance,
          threshold = threshold[ivar]
          # Number.Edges.per.Point = npositives
        )
      )

    } #end-function-constructor



  H <- constructHOMOLOGY(ivar, Ydat, Xdat, dimension, threshold, method, angle)
      mp_union <- H[["manifold_unioned"]]

      mp_reflection <- estimate_symmetric_reflection(mp_union)

      bb <- sf::st_make_grid(x = mp_union, n = 1)

      mp_sym_difference <-
        sf::st_sym_difference(mp_union, mp_reflection)


      H$sym_difference <- mp_sym_difference
      H$sym_reflection <- mp_reflection

      Symmetric.Diff.Area <- sf::st_area(mp_sym_difference)
      Manifold.Area <- sf::st_area(mp_union)
      Box.Area <- sf::st_area(bb)

      return(
        list(
          threshold = H[["threshold"]],
          Number.Edges.per.Point = H[["Number.Edges.per.Point"]],
          Neigborhood.Distance = H[["neigborhood.distance"]],
          Manifold.Area = Manifold.Area,
          Box.Area = Box.Area,
          Geometric.R2 = 1 - Manifold.Area / Box.Area,
          Symmetric.Diff.Area = Symmetric.Diff.Area,
          Geometric.Sensitivity = Symmetric.Diff.Area / (2 * Manifold.Area),
          homology = H
        )
      )


}

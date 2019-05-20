#' Title
#'
#' @param Ydat numeric vector of responses in a model.
#' @param Xdat numeric matrix or data.frame of covariables.
#' @param dimension number of homology groups to estimate. The only value
#'   accepted is \code{dimension = 3}, but in future release this will change.
#' @param threshold percent of radius or sizes of triangles to keep to smooth
#'   the homology complex. Defaults to `3`.
#' @param knearest number of the nearest neigbors keep to built the homology
#'   complex
#' @param method type of method to build the homology complex. Two choices are
#'   accepted: \code{delanauy} o \code{VR} (Vietoris-Rips).
#' @param mc.cores number of cores to use in the procedure. Defaults \code{1}.
#'
#' @return Put output
#'
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
#' estimation <- TopSA(Ydat = Y, Xdat = X)


#'
#' @importFrom methods is
#' @importFrom stats dist quantile
TopSA <-
  function(Ydat,
           Xdat,
           dimension = 3,
           threshold = 0.05,
           knearest = 20,
           method = c("delanauy", "VR"),
           mc.cores = 1) {
    #Arguments:
    #Y: matrix of model outputs (only one column)
    #X: matrix model parameters
    #radius: radius to build the neighborhood graph
    #dimension: number of homology spaces
    #MAIN
    if (!requireNamespace("scales", quietly = TRUE)) {
      stop("Please install the package scales: install.packages('scales')")
    }#end-require-scales
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop("Please install the package igraph: install.packages('igraph')")
    }#end-require-igraph
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("Please install the package sf: install.packages('sf')")
    }#end-require-sf
    # # if (!requireNamespace("rgeos", quietly = TRUE)) {
    # #   stop("Please install the package rgeos: install.packages('rgeos')")
    # # }#end-require-rgeos
    # if (!requireNamespace("pbmcapply", quietly = TRUE)) {
    #   stop("Please install the package pbmcapply: install.packages('pbmcapply')")
    # }#end-require-pbmcapply
    # # if (!requireNamespace("multimode", quietly = TRUE)) {
    # #   stop("Please install the package multimode: install.packages('multimode')")
    # # }#end-require-multimode


    Xdat <- as.data.frame(Xdat)
    Ydat <- as.data.frame(Ydat)

    ANS <- list()
    ANS[['call']] <- match.call()
    ANS[['Xdat']] <- Xdat
    ANS[['Ydat']] <- Ydat
    ANS[['dimension']] <- dimension
    ANS[['threshold']] <- threshold
    par.names = colnames(Xdat)	#gets parameters names
    if (is.null(colnames(Xdat)))
      par.names = paste0('X', 1:ncol(Xdat))

    message("Index estimation")


    sensitivity_results <- try(pbmcapply::pbmclapply(
      X = 1:ncol(Xdat),
      FUN = estimate_sensitivity_index,
      Ydat = Ydat,
      Xdat = Xdat,
      dimension = dimension,
      knearest = knearest,
      threshold = threshold,
      method = method,
      mc.cores.inner = mc.cores,
      mc.cores = mc.cores
    ),
    silent = TRUE)


    for (i in 1:ncol(Xdat)) {
      sensitivity_results[[i]]$xname <- colnames(Xdat)[i]
      sensitivity_results[[i]]$yname <- colnames(Ydat)
    }

    ANS[["results"]] <- sensitivity_results
    class(ANS) <- 'TopSA'
    return(ANS)

    # SA.tab <- try(pbmcapply::pbmclapply(X = HOMOLOGY,
    #                                     FUN = estimate_index,
    #                                     mc.cores = nc),
    #               silent = T
    # )
    #
    # if (is(SA.tab, 'try-error')) {
    #Windows does not support mclapply...

    # }

    #  sensitivity_table <- t(sapply(sensitivity_result, function(x) {
    #   as.numeric(x[1:5])
    # }))
    #
    # colnames(sensitivity_table) <- c('Radius', 'Manifold Area', 'Box Area' , 'Geometric correlation', 'Symmetric index')
    # rownames(sensitivity_table) <- par.names
    # ANS[['sensitivity_table']] <- sensitivity_table

  }

estimate_sensitivity_index <-
  function (ivar,
            Ydat,
            Xdat,
            dimension,
            knearest,
            threshold,
            method,
            mc.cores.inner) {
    if (method == "delanauy") {
      index_Obj <-
        Delanauy_homology(ivar, Ydat, Xdat, dimension, threshold)
    } else if (method == "VR") {
      index_Obj <-  VR_homology(ivar, Ydat, Xdat, dimension, knearest, threshold, mc.cores.inner)
    } else{
      index_Obj <- NULL
      stop("No method defined")
    }
    return(index_Obj)
  }




VR_homology <- function(ivar, Ydat, Xdat, dimension, knearest, threshold, mc.cores.inner) {
  constructHOMOLOGY <-
    function (ivar, Ydat, Xdat, dimension, threshold) {
      Y <- as.matrix(Ydat)
      X <- as.matrix(Xdat[, ivar])

      Xr <- scales::rescale(X , to = c(0, 1))
      Yr <- scales::rescale(Y , to = c(0, 1))
      neigborhood.distance <- dist(cbind(Xr, Yr))

      #radius[i] <- fivenum(neigborhood.distance)[2]
      # dY <- diff(range(Y))
      # dX <- diff(range(X))
      nd <- as.numeric(neigborhood.distance)

      r <- quantile(nd, threshold)

      #r <- as.numeric(fivenum(neigborhood.distance)[2])
      #  if (dY >= dX) {
      #    s <- dX / dY
      #    r <-  m * s
      # } else{
      #   s <- dY / dX
      #    r <- m * s
      #  }


      adjacency.matrix <-
        as.matrix(neigborhood.distance) * (as.matrix(neigborhood.distance) <= r)
      adjacency.matrix[adjacency.matrix==0] <- NA

      npositives <- NULL
      for (k in 1:length(Xr)) {
        npositives[k] <- sum(adjacency.matrix[k,] >0,na.rm = TRUE)
        knearestidx <- order(adjacency.matrix[k,])[1:knearest]
        idx <- 1:length(Xr) %in% knearestidx
        adjacency.matrix[k, !idx] <- NA
      }

      adjacency.matrix[is.na(adjacency.matrix)]<- 0
      adjacency.matrix <- adjacency.matrix > 0




      graphBase <-
        igraph::graph.adjacency(adjacency.matrix,
                                mode = "undirected",
                                weighted = NULL)
      graphBase <-
        igraph::set_vertex_attr(graphBase, name = "x", value = X)
      graphBase <-
        igraph::set_vertex_attr(graphBase, name = "y", value = Y)

      H <- list()

      H[[1]] <-
        as.integer(igraph::as_data_frame(graphBase, "vertices")[, "name"])
      v <-
        igraph::as_data_frame(graphBase, "edges")[, c("to", "from")]
      v <-  data.matrix(v)
      mode(v) <- "integer"
      H[[2]] <-  v

      for (d in 3:dimension) {
        clq <- igraph::cliques(graphBase, min = d, max = d)
        H[[d]] <- do.call("rbind", clq)
      }

      return(
        list(
          graph = graphBase,
          homology = H,
          neigborhood.distance = neigborhood.distance,
          threshold = r,
          Number.Edges.per.Point = npositives
        )
      )


    } #end-function-constructor

  Ydat <- as.matrix(Ydat)

  H <- constructHOMOLOGY(ivar, Ydat, Xdat, dimension, threshold)

  vv <- igraph::as_data_frame(H[["graph"]], "vertices")
  H2 <- H[["homology"]][[3]]
  l <- list()
  idxObj <- sort(unique(as.numeric(H2)))
  Vertices <- data.frame (name = vv$name[idxObj],
                          x = vv$x[idxObj],
                          y = vv$y[idxObj])

  l <-  try(pbmcapply::pbmclapply(
    X = 1:nrow(H2),
    FUN = function(i) {
      idxTriangle <- H2[i, ]
      idxname <- Vertices$name %in% c(idxTriangle, idxTriangle[1])
      subVertices <-  Vertices[idxname ,]
      subVertices <- rbind(subVertices, subVertices[1, ])
      Triangle <- as.matrix((subVertices[, c(2, 3)]))
      rownames(Triangle) <- NULL
      p <- sf::st_polygon(list(Triangle))
    },
    mc.cores = mc.cores.inner
  ),
  silent = T)

  if (is(l, 'try-error')) {
    #Windows does not support mclapply...
    l <- lapply(
      X = 1:nrow(H2),
      FUN = function(i) {
        idxTriangle <- H2[i, ]
        idxname <- Vertices$name %in% c(idxTriangle, idxTriangle[1])
        subVertices <-  Vertices[idxname ,]
        subVertices <- rbind(subVertices, subVertices[1, ])
        Triangle <- as.matrix((subVertices[, c(2, 3)]))
        rownames(Triangle) <- NULL
        p <- sf::st_polygon(list(Triangle))
      }
    )
  }

  mp <- sf::st_multipolygon(l)
  mp_union <- sf::st_union(mp) - c(min(Xdat[, ivar]), min(Ydat))
  bb <- sf::st_make_grid(x = mp_union, n = 1)

  reflectiony <-  matrix(c(1, 0, 0, -1), 2, 2)

  mp_reflectiony <-
    mp_union * reflectiony + c(0,   min(mp_union[[1]][, 2]) +  max(mp_union[[1]][, 2]))


  # plot(mp_union, asp = 0, col = "blue", axes = TRUE)
  # plot(mp_reflectiony,
  #      asp = 0,
  #      col = "red",
  #      add = TRUE)


  mp_sym_difference <-
    sf::st_sym_difference(mp_union, mp_reflectiony)



  # plot(mp_sym_difference, col = "blue", asp=0)
  mp_sym_difference_area <- sf::st_area(mp_sym_difference)

  Manifold.Area <- sf::st_area(mp_union)
  Box.Area <- sf::st_area(bb)

  return(
    list(
      threshold = H[["threshold"]],
      Number.Edges.per.Point = H[["Number.Edges.per.Point"]],
      Manifold.Area = Manifold.Area,
      Box.Area = Box.Area,
      Geometric.Correlation = 1 - Manifold.Area / Box.Area,
      Symmetric.Diff.Area = mp_sym_difference_area,
      Symmetric.Index = mp_sym_difference_area / (2 * Manifold.Area),
      manifold.plot = mp_union + c(min(Xdat[, ivar]), min(Ydat))
    )
  )
}

Delanauy_homology <-
  function(ivar, Ydat, Xdat, dimension, threshold) {

    idx <- order(Xdat[,ivar])
    Xdat[,ivar] <- Xdat[idx,ivar]
    Ydat <- as.matrix(Ydat[idx,])

    rx <- range(Xdat[,ivar])
    ry <- range(Ydat)
    minx <- min(Xdat[,ivar])
    miny <- min(Ydat)
    maxy <- max(Ydat)
    datapoints <-
      sf::st_as_sf(x = data.frame(
        x = (Xdat[, ivar] - minx) / (diff(rx)),
        y = (Ydat - miny) / diff(ry)
      ),
      coords = c("x", "y"))

    triangulation <-
      sfdct::ct_triangulate(sf::st_union(datapoints),
                            S = length(datapoints$geometry),
                            q = 40)
    single_triangles <- sf::st_collection_extract(triangulation)
    # plot(
    #   single_triangles,
    #   axes = TRUE,
    #   asp = 1,
    #   col = "blue",
    #   border = "white"
    # )
    single_triangles_area <-  sf::st_area(single_triangles)
    # plot(single_triangles_area)
    # abline(h = boxplot(single_triangles_area, plot = FALSE)$stats,
    #        col = "red")
    # abline(h = mean(single_triangles_area) + sd(single_triangles_area),
    #        col = "blue")
    # abline(h = median(single_triangles_area) , col = "green")
    if (exists("threshold") == FALSE) {
     threshold <- 0.90
    }
    cutoff <- quantile(single_triangles_area, threshold)
    idx <- single_triangles_area <=cutoff

    single_triangles <- sf::st_multipolygon(single_triangles[idx])

    A <- matrix(c(diff(rx),0,0,diff(ry)),nrow = 2)
    #b <- c(minx, miny)
    single_triangles <- single_triangles*A #+b
    # plot(single_triangles,
    #      asp = 1,
    #      col = "blue",
    #      border = "white")

    mp_union <-
      sf::st_union(single_triangles)

    bb <- sf::st_make_grid(x = mp_union, n = 1)

    bbox <- sf::st_bbox(mp_union)

    reflectiony <-  matrix(c(1, 0, 0, -1), 2, 2)

    pointsInComplex <- do.call(rbind, do.call(rbind, mp_union))

    mp_reflectiony <-
      (mp_union)  * reflectiony +
      c(0,   min(pointsInComplex[, 2]) + max(pointsInComplex[, 2]))

    mp_sym_difference <-
      sf::st_sym_difference(mp_union, mp_reflectiony)



    mp_sym_difference_area <- sf::st_area(mp_sym_difference)

    Manifold.Area <- sf::st_area(mp_union)
    Box.Area <- sf::st_area(bb)

    return(
      list(
        threshold = cutoff,
        Manifold.Area = Manifold.Area,
        Box.Area = Box.Area,
        Geometric.Correlation = 1 - Manifold.Area / Box.Area,
        Symmetric.Index = mp_sym_difference_area / (2 * Manifold.Area),
        manifold.plot = mp_union + c(min(Xdat[, ivar]), min(Ydat))
      )
    )
  }







  # ggplot2::ggplot() +
  #   ggplot2::geom_sf(data = manifold, ggplot2::aes(fill = "Estimated")) +
  #   ggplot2::geom_sf(data = manifold_reflectionx, ggplot2::aes(fill = "Sym. Reflection")) +
  #   ggplot2::geom_sf(data = manifold_intersection, ggplot2::aes(fill = "Intersection")) +
  #   ggplot2::scale_fill_manual(
  #     name = "Legend",
  #     values = c(
  #       "Estimated" = "deepskyblue",
  #       "Sym. Reflection" =  "orange",
  #       "Intersection" = "grey90"
  #     ),
  #     breaks = c("Estimated" ,
  #                "Sym. Reflection",
  #                "Intersection")
  #   ) +
  #   ggplot2::geom_sf(data = datapoints, shape = "*", size = 2) +
  #   ggplot2::geom_sf(
  #     data = boundingbox,
  #     color = "red",
  #     fill = NA,
  #     size = 1
  #   ) +
  #   ggplot2::theme_minimal()





  #   # plotdata <-
  # #   data.frame(
  # #     TypeObject = c("Original", "Reflection", "Sym Diff"),
  # #     AreaObject = c()
  # #     geom = c(manifold, manifold_reflectionx, manifold_sym_difference)
  # #   )
  #
  #
  # if (symmetric.diff) {
  #   ggplot2::ggplot() +
  #     ggplot2::geom_sf(data = manifold_sym_difference, fill = "orange") +
  #     ggplot2::geom_sf(data = datapoints,
  #             shape = "*",
  #             size = 2) +
  #     ggplot2::geom_sf(
  #       data = boundingbox,
  #       color = "red",
  #       fill = NA,
  #       size = 1
  #     ) +
  #     ggplot2::theme_minimal()
  # } else{
  #   ggplot2::ggplot() +
  #     ggplot2::geom_sf(data = manifold, fill = "deepskyblue", color = "blue" ) +
  #     ggplot2::geom_sf(data = datapoints, shape = "*", size = 2) +
  #     ggplot2::geom_sf(data = boundingbox, color = "red", fill = NA, size = 1 ) +
  #     ggplot2::theme_minimal()
  #
  #   if (with.reflection == FALSE & legend == TRUE) {
  #
  #     #     legend(
  #     #   x = "bottomright",
  #     #   y.intersp = 1,
  #     #   legend = c(
  #     #     paste0(
  #     #       "Manifold Area: ",
  #     #       round(TopSAObj$results[[nvar]]$Manifold.Area, 2)
  #     #     ),
  #     #     paste0("Box Area: ", round(TopSAObj$results[[nvar]]$Box.Area, 2)),
  #     #     paste0(
  #     #       "Geometric Correlation: ",
  #     #       round(TopSAObj$results[[nvar]]$Geometric.Correlation, 2)
  #     #     )
  #     #   ),
  #     #   border = c("blue", "red", "white"),
  #     #   fill = c("deepskyblue", "white", "white")
  #     # )
  #
  #   }
  #
  #   if (with.reflection == TRUE) {
  #
  #     ggplot2::ggplot() +
  #       ggplot2::geom_sf(data = manifold, fill = "deepskyblue" ) +
  #       ggplot2::geom_sf(data = manifold_reflectionx, fill = "seagreen1" ) +
  #       ggplot2::geom_sf(data = datapoints, shape = "*", size = 2) +
  #       ggplot2::geom_sf(data = boundingbox, color = "red", fill = NA, size = 1 ) +
  #       ggplot2::theme_minimal()
  #
  #
  #
  #     # plot(manifold_reflectionx, col = "seagreen1", add = TRUE)
  #     # plot(
  #     #   datapoints,
  #     #   col = "black",
  #     #   type = "p",
  #     #   pch = ".",
  #     #   cex = 3,
  #     #   add = TRUE
  #     # )
  #
  #
  #     if (legend == TRUE) {
  #       legend(
  #         x = "bottomright",
  #         y.intersp = 1,
  #         legend = c(
  #           paste0(
  #             "Manifold Area: ",
  #             round(TopSAObj$results[[nvar]]$Manifold.Area, 2)
  #           ),
  #           paste0(
  #             "Symmetric Reflection Area: ",
  #             round(TopSAObj$results[[nvar]]$Manifold.Area, 2)
  #           ),
  #           paste0("Box Area: ", round(TopSAObj$results[[nvar]]$Box.Area, 2)),
  #           paste0(
  #             "Geometric Correlation: ",
  #             round(TopSAObj$results[[nvar]]$Geometric.Correlation, 2)
  #           )
  #         ),
  #         border = c("blue", "black", "red", "white"),
  #         fill = c("deepskyblue", "seagreen1", "white", "white")
  #       )
  #     }
  #   }
  # }
#}


# plot.TopSA <- function(HOMOLOGYObj, n = 5000,...) {
#   if (!requireNamespace("igraph", quietly = TRUE)) {
#     stop("Please install the package np: install.packages('igraph')")
#   }#end-require-igraph
#   g <- HOMOLOGYObj$graph
#   H2 <- HOMOLOGYObj$homology[[3]]
#   bb <- HOMOLOGYObj$box
#
#   par(cex = 2)
#   plot.igraph(
#     g,
#     rescale = FALSE,
#     axes = TRUE,
#     asp = 0,
#     ylim = range(V(g)$y[as.integer(H2)]),
#     xlim = range(V(g)$x[as.integer(H2)]),
#     vertex.size = 0,
#     vertex.label = NA,
#     vertex.label.cex = 1.5,
#     edge.arrow.size = 1,
#     edge.width = 0.5,
#     edge.color = "lightgrey",
#     ...
#   )
#
#
#   if (n == "full") {
#     s <- 1:nrow(H2)
#   } else{
#     s <- sample(1:nrow(H2), min(nrow(H2), n))
#   }
#
#   H2 <- H2[s,]
#   for (k in 1:dim(H2)[1]) {
#     polygon(
#       V(g)$x[H2[k,]],
#       V(g)$y[H2[k,]],
#       border = adjustcolor("blue", alpha.f = 0.1),
#       col = adjustcolor("lightblue", alpha.f = 0.9)
#     )
#   }
#
#   rect(
#     xleft = bb[1, 1],
#     xright = bb[1, 2],
#     ybottom = bb[2, 1],
#     ytop = bb[2, 2],
#     lwd = 2,
#     border = adjustcolor("red", alpha.f = 0.8)
#   )
#
#
# }

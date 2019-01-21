# Author: Maikol Solís
#To create a robust function, I took as a template the script of Filippo Monari,
#'sobolSmthSpl' from the package 'sensitivity'. However, the core estimations
#are my original work.

TopSA <-
  function(Ydat,
           Xdat,
           dimension = 3,
           threshold = 0.05,
           method = c("delanauy", "VR")) {
    #Arguments:
    #Y: matrix of model outputs (only one column)
    #X: matrix model parameters
    #radius: radius to build the neighborhood graph
    #dimension: number of homology spaces
    #MAIN
    # if (!requireNamespace("scales", quietly = TRUE)) {
    #   stop("Please install the package scales: install.packages('scales')")
    # }#end-require-scales
    # if (!requireNamespace("igraph", quietly = TRUE)) {
    #   stop("Please install the package igraph: install.packages('igraph')")
    # }#end-require-igraph
    # if (!requireNamespace("sf", quietly = TRUE)) {
    #   stop("Please install the package sf: install.packages('sf')")
    # }#end-require-sp
    # # if (!requireNamespace("rgeos", quietly = TRUE)) {
    # #   stop("Please install the package rgeos: install.packages('rgeos')")
    # # }#end-require-rgeos
    # if (!requireNamespace("pbmcapply", quietly = TRUE)) {
    #   stop("Please install the package pbmcapply: install.packages('pbmcapply')")
    # }#end-require-pbmcapply
    # # if (!requireNamespace("multimode", quietly = TRUE)) {
    # #   stop("Please install the package multimode: install.packages('multimode')")
    # # }#end-require-multimode



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
      threshold = threshold,
      method = method,
      mc.cores = min(ncol(Xdat), parallel::detectCores() - 1)
    ),
    silent = TRUE

    )


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
            threshold,
            method) {
    if (method == "delanauy") {

    } else if (method == "VR") {
      index_Obj <-  VR_homology(ivar, Ydat, Xdat, dimension, threshold)
    } else{
      stop("No method defined")
    }
    return(index_Obj)
  }




VR_homology <- function(ivar, Ydat, Xdat, dimension, threshold) {

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


      adjacency.matrix <- as.matrix(neigborhood.distance) <= r
      diag(adjacency.matrix) <- 0

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
          radius = r
        )
      )


    } #end-function-constructor


  H <- constructHOMOLOGY(ivar, Ydat, Xdat, dimension, threshold)

  vv <- igraph::as_data_frame(H[["graph"]], "vertices")
  H2 <- H[["homology"]][[3]]
  l <- list()
  idxObj <- sort(unique(as.numeric(H2)))
  Vertices <- data.frame (name = vv$name[idxObj],
                          x = vv$x[idxObj],
                          y = vv$y[idxObj])

  nc <-  ceiling((parallel::detectCores() - 1 - ncol(Xdat)) / ncol(Xdat))

  l <-  try(pbmcapply::pbmclapply(
    X = 1:nrow(H2),
    FUN = function(i) {
      idxTriangle <- H2[i,]
      idxname <- Vertices$name %in% c(idxTriangle, idxTriangle[1])
      subVertices <-  Vertices[idxname , ]
      subVertices <- rbind(subVertices, subVertices[1,])
      Triangle <- as.matrix((subVertices[, c(2, 3)]))
      rownames(Triangle) <- NULL
      p <- sf::st_polygon(list(Triangle))
    },
    mc.cores = nc
  ),
  silent = T)

  if (is(l, 'try-error')) {
    #Windows does not support mclapply...
    l <- lapply(
      X = 1:nrow(H2),
      FUN = function(i) {
        idxTriangle <- H2[i,]
        idxname <- Vertices$name %in% c(idxTriangle, idxTriangle[1])
        subVertices <-  Vertices[idxname , ]
        subVertices <- rbind(subVertices, subVertices[1,])
        Triangle <- as.matrix((subVertices[, c(2, 3)]))
        rownames(Triangle) <- NULL
        p <- sf::st_polygon(list(Triangle))
      }
    )
  }

  mp <- sf::st_multipolygon(l)
  mp_union <- sf::st_union(mp)
  bb <- sf::st_make_grid(x = mp_union, n = 1)

  reflectiony <-  matrix(c(1, 0, 0, -1), 2, 2)

  mp_reflectiony <-
    mp_union * reflectiony + c(0, 2 * sf::st_centroid(mp)[2])

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
      radius = H[["radius"]],
      Manifold.Area = Manifold.Area,
      Box.Area = Box.Area,
      Geometric.Correlation = 1 - Manifold.Area / Box.Area,
      Symmetric.Index = mp_sym_difference_area / (2 * Manifold.Area),
      manifold.plot = mp
    )
  )
}








print.TopSA <- function(x, ...) {

   sensitivity_table <- t(sapply(x$results, function(x) {
    as.numeric(x[1:5])
  }))

  colnames(sensitivity_table) <- c('Radius', 'Manifold Area', 'Box Area' , 'Geometric correlation', 'Symmetric index')
  rownames(sensitivity_table) <- par.names


  cat("\nCall:\n", deparse(x[['call']]), "\n", sep = "")
  cat("\nNumber of variables:", ncol(x[['Xdat']]), "\n")
  cat("\nNumber of observations:", length(x[['Ydat']]), "\n")
  cat("\nFirst order indices\n")
  print(sensitivity_table[,])
  # return(sensitivity_table)
}



plot.TopSA <- function(resultsObj, percent = 0.1, ...) {
  if (percent < 0 | percent > 1)
    stop("Percent is between 0 and 1")

  sapply(resultsObj$results, function(p) {
    p <- p$manifold.plot
    bb <- sf::st_make_grid(x = p, n = 1)
    idx <- sample(1:length(p), length(p) * percent)
    p <- sf::st_multipolygon(p[idx])

    plot(
      p,
      asp = 0,
      axes = TRUE,
      col = "blue",
      border = "lightblue"
    )
    plot(bb, add = TRUE)
  })
  invisible(NULL)
}




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

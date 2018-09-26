# Author: Maikol Solís
#To create a robust function, I took as a template the script of Filippo Monari,
#'sobolSmthSpl' from the package 'sensitivity'. However, the core estimations
#are my original work.

sobolManifold <-
  function(Y,
           X,
           radius = 1,
           dimension = 3,
           alpha = 0.2) {
    #Arguments:
    #Y: matrix of model outputs (only one column)
    #X: matrix model parameters
    #radius: radius to build the neighborhood graph
    #dimension: number of homology spaces
    #MAIN

    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop("Please install the package np: install.packages('igraph')")
    }#end-require-igraph
    if (!requireNamespace("sp", quietly = TRUE)) {
      stop("Please install the package sp: install.packages('sp')")
    }#end-require-sp
    if (!requireNamespace("rgeos", quietly = TRUE)) {
      stop("Please install the package rgeos: install.packages('rgeos')")
    }#end-require-rgeos
    if (!requireNamespace("pbmcapply", quietly = TRUE)) {
      stop("Please install the package pbmcapply: install.packages('pbmcapply')")
    }#end-require-pbmcapply
    # if (!requireNamespace("multimode", quietly = TRUE)) {
    #   stop("Please install the package multimode: install.packages('multimode')")
    # }#end-require-multimode



    ANS <- list()
    ANS[['call']] <- match.call()
    ANS[['X']] <- X
    ANS[['Y']] <- Y
    ANS[['dimension']] <- dimension
    ANS[['alpha']] <- alpha
    par.names = colnames(X)	#gets parameters names
    if (is.null(colnames(X)))
      par.names = paste0('X', 1:ncol(X))
    message("Homology construction")
    HOMOLOGY <- constructHOMOLOGY(Y, X, radius, dimension, alpha)

    # nc <-  min(ncol(X), parallel::detectCores()) / 2
    #
    # SA.tab <- try(pbmcapply::pbmclapply(X = HOMOLOGY,
    #                                     FUN = estimate_index,
    #                                     mc.cores = nc),
    #               silent = T
    # )
    #
    # if (is(SA.tab, 'try-error')) {
    #   #Windows does not support mclapply...
    message("Index estimation")
    SA.tab <- lapply(X = HOMOLOGY,
                     FUN = estimate_index)
    # }

    axis.scale <- NULL
    for (i in 1:ncol(X)) {
      HOMOLOGY[[i]]$box <- SA.tab[[i]]$box
      radius[i] <- HOMOLOGY[[i]]$radius
      # axis.scale[i] <- HOMOLOGY[[i]]$scales
    }
    ANS[['radius']] <- radius
    # ANS[['axis.scale']] <- axis.scale

    ANS[['HOMOLOGY']] <- HOMOLOGY
    SA.print <- t(sapply(SA.tab, function(x) {
      as.numeric(x[1:3])
    }))
    colnames(SA.print) <- c('Obj Area', 'Square Area' , 'Index')
    rownames(SA.print) <- par.names
    ANS[['index']] <- SA.print
    class(ANS) <- 'sobolManifold'
    return(ANS)
  }

estimate_index <- function (H) {
  vv <- igraph::as_data_frame(H[["graph"]], "vertices")
  H2 <- H[["homology"]][[3]]
  l <- list()
  idxObj <- sort(unique(as.numeric(H2)))
  Vertices <- data.frame (name = vv$name[idxObj],
                          x = vv$x[idxObj],
                          y = vv$y[idxObj])

  nc <-   parallel::detectCores() - 1

  l <-  try(pbmcapply::pbmclapply(
    X = 1:nrow(H2),
    FUN = function(i) {
      idxTriangle <- H2[i,]
      idxname <- Vertices$name %in% c(idxTriangle, idxTriangle[1])
      subVertices <-  Vertices[idxname , ]
      subVertices <- rbind(subVertices, subVertices[1,])
      Triangle <- as.matrix((subVertices[, c(2, 3)]))
      p <- sp::Polygon(Triangle, hole = FALSE)
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
        p <- sp::Polygon(Triangle, hole = FALSE)
      }
    )
  }

  ps <- sp::Polygons(l, 1)
  sps <- sp::SpatialPolygons(list(ps))
  sps <- rgeos::gUnion(sps, sps)

  ObjArea <- rgeos::gArea(sps)
  bb <- sp::bbox(sps)
  SqArea <- prod(diff(t(bb)))

  return(list(
    ObjArea = ObjArea,
    SqArea = SqArea,
    index = 1 - ObjArea / SqArea,
    box = bb
  ))

}

constructHOMOLOGY <- function (Y, X, radius, dimension, alpha) {
  library(igraph)
  #DOC
  #Arguments:
  #Y: matrix containing the response variable
  #X: matrix of inputs
  #radius: radius to build the neighborhood graph
  #dimension: number of homology spaces

  #CONTAINS

  lower_nbrs <- function(graph, node) {
    n <- igraph::neighbors(graph, node)
    edgesreturn <-  as.numeric(n[n < node])
    # nodeSet =  as.numeric(V(graph))
    # edgesreturn <- NULL
    # for (x in nodeSet) {
    #   if (length(E(induced_subgraph(graph, c(x, node)))) == 1 &
    #       node > x) {
    #     edgesreturn <- c(edgesreturn, x)
    #   }
    # }
    return(edgesreturn)
  }

  add_cofaces <- function(graph, k, tau, N, V) {
    V <- union(V, list(tau))
    if (length(tau) >= k) {
      return(V)
    } else{
      for (v in N) {
        sigma <- union(tau, v)
        M <- intersect(N, lower_nbrs(graph = graph, v))
        V <- add_cofaces(graph, k, sigma, M, V)
      }
    }
    return(V)
  }

  incremental_VR <- function(graph, k) {
    nodeSet =  as.numeric(V(graph))
    V <- list()
    for (u in nodeSet) {
      N <- lower_nbrs(graph, u)
      V <- add_cofaces(graph, k, u, N, V)
    }
    return(V)
  }

  extract_vertices <- function(x, k) {
    if (!is.null(x) & length(x) == k) {
      return(x)
    }
  }

  create_homology_groups <- function(k, VR) {
    H <- sapply(VR, extract_vertices , k)
    H <- H[!sapply(H, is.null)]
    H <- do.call('rbind', H)
  }

  constructor <-
    function(i, Ydat, Xdat, radius, dimension, alpha) {
      Y <- Ydat
      X <- Xdat

      Xr <- scales::rescale(X , to = c(0, 1))
      Yr <- scales::rescale(Y , to = c(0, 1))
      neigborhood.distance <- dist(cbind(Xr[, i], Yr))

      #radius[i] <- fivenum(neigborhood.distance)[2]
      # dY <- diff(range(Y))
      # dX <- diff(range(X))
      nd <- as.numeric(neigborhood.distance)

      r <- quantile(nd, alpha)

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
        igraph::set_vertex_attr(graphBase, name = "x", value = X[, i])
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


      # VR <- incremental_VR(graph = graphBase, k = dimension)

      # H <- lapply(X = 1:dimension, FUN = create_homology_groups, VR)


      return(
        list(
          graph = graphBase,
          homology = H,
          neigborhood.distance = neigborhood.distance,
          radius = r
        )
      )


    } #end-function-constructor

  nc <-  min(ncol(X), parallel::detectCores())

  ANS <- try(pbmcapply::pbmclapply(
    X = 1:ncol(X),
    FUN = constructor,
    Ydat = Y,
    Xdat = X,
    radius = radius,
    dimension = dimension,
    alpha = alpha,
    mc.cores = nc
  ),
  silent = T)

  if (is(ANS, 'try-error')) {
    #Windows does not support mclapply...
    HOMOLOGY <- lapply(
      X = 1:ncol(X),
      FUN = constructor,
      Ydat = Y,
      Xdat = X,
      radius = radius,
      dimension = dimension,
      alpha = alpha
    )
  } else {
    HOMOLOGY <<- ANS
  }

  return(HOMOLOGY)
} #end-optNP


print.sobolManifold <- function(x, ...) {
  cat("\nCall:\n", deparse(x[['call']]), "\n", sep = "")
  cat("\nNumber of observations:", length(x[['Y']]), "\n")
  cat("\nFirst order indices\n")
  print(x[['index']])
}


plot.sobolManifold <- function(HOMOLOGYObj, n = 5000) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Please install the package np: install.packages('igraph')")
  }#end-require-igraph
  g <- HOMOLOGYObj$graph
  H2 <- HOMOLOGYObj$homology[[3]]
  bb <- HOMOLOGYObj$box

  par(cex = 2)
  plot.igraph(
    g,
    rescale = FALSE,
    axes = TRUE,
    asp = 0,
    ylim = range(V(g)$y[as.integer(H2)]),
    xlim = range(V(g)$x[as.integer(H2)]),
    vertex.size = 0,
    vertex.label = NA,
    vertex.label.cex = 1.5,
    edge.arrow.size = 1,
    edge.width = 0.5,
    edge.color = "lightgrey"
  )


  if (n == "full") {
    s <- 1:nrow(H2)
  } else{
    s <- sample(1:nrow(H2), min(nrow(H2), n))
  }

  H2 <- H2[s,]
  for (k in 1:dim(H2)[1]) {
    polygon(
      V(g)$x[H2[k,]],
      V(g)$y[H2[k,]],
      border = adjustcolor("blue", alpha.f = 0.1),
      col = adjustcolor("lightblue", alpha.f = 0.9)
    )
  }

  rect(
    xleft = bb[1, 1],
    xright = bb[1, 2],
    ybottom = bb[2, 1],
    ytop = bb[2, 2],
    lwd = 2,
    border = adjustcolor("red", alpha.f = 0.8)
  )


}

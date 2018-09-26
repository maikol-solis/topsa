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
  H2 <- HOMOLOGY[["homology"]][[3]]
  idxObj <- sort(unique(as.numeric(H2)))
  Vertices <- data.frame (name = vv$name[idxObj],
                          x = normalizeX(vv$x[idxObj]),
                          y = normalizeX(vv$y[idxObj]))

  for (i in 1:nrow(H2)) {
    idxTriangle <- H2[i,]
    idxname <- Vertices$name %in% c(idxTriangle, idxTriangle[1])
    subVertices <-  Vertices[idxname , ]
    subVertices <- rbind(subVertices, subVertices[1,])
    Triangle <- as.matrix((subVertices[, c(2, 3)]))
    p <- sp::Polygon(Triangle, hole = FALSE)
    l[[i]] <- p
  }
  ps <- sp::Polygons(l,1)
  sps <- sp::SpatialPolygons(list(ps))
  sps <- rgeos::gUnion(sps, sps)

  AreaObj <- rgeos::gArea(sps)
  bb <- sp::bbox(sps)
  AreaCuad <- prod(diff(t(bb)))

  return(c(AreaObj, AreaCuad, 1 - AreaObj / AreaCuad))

  #calculates the standard error of the main effect estimates
  # yi.sc = (yi - mean(yi)) / sd(yi)            	#scaled yi
  # u = (yi - gi) / (sd(yi) * abs(1 - Si) ** 0.5)   #scales residuals
  # Si.se = abs(1 - Si) * sd(yi.sc ** 2 - u ** 2) / length(yi) ** 0.5
  # q0.05 = qnorm(0.05, Si, Si.se)
  # return(c(Si, bw, Si.se, q0.05))
}

constructHOMOLOGY <- function (Y, X, radius, dimension) {
  library(igraph)
  #DOC
  #Estimates nonparametrically the regression curve across all the column of the
  #matrix Y.
  #Arguments:
  #Y: matrix containing the response variable
  #X: matrix of inputs
  #radius: radius to build the neighborhood graph
  #dimension: number of homology spaces

  #CONTAINS

  lower_nbrs <- function(graph, node) {
    nodeSet =  as.numeric(V(graph))
    edgesreturn <- NULL
    for (x in nodeSet) {
      if (length(E(induced_subgraph(graph, c(x, node)))) == 1 &
          node > x) {
        edgesreturn <- c(edgesreturn, x)
      }
    }
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

  create_homology_groups <- function(k,VR) {
    H <- sapply(VR, extract_vertices , k)
    H <- H[!sapply(H, is.null)]
    H <- do.call('rbind', H)
  }

  constructor <-
    function(i, Ydat,Xdat,radius,dimension) {
      Y <- Ydat
      X <- Xdat

      neigborhood.distance <- as.matrix(dist(cbind(X[, i], Y[, i])))
      adjacency.matrix <-
        neigborhood.distance * (neigborhood.distance < radius)
      diag(adjacency.matrix) <- 0

      graphBase <-
        igraph::graph.adjacency(adjacency.matrix,
                                mode = "undirected",
                                weighted = TRUE)
      graphBase <-
        igraph::set_vertex_attr(graphBase, name = "x", value = X[, i])
      graphBase <-
        igraph::set_vertex_attr(graphBase, name = "y", value = Y[, i])

      VR <- incremental_VR(graph = graphBase, k = dimension)

      H <- lapply(X = 1:dimension, FUN = create_homology_groups, VR)


      return(list(graph = graphBase, homology = H))


    } #end-function-constructor

  nc <-  min(ncol(Y), parallel::detectCores())

  ANS <- try(parallel::mclapply(
    X = 1:ncol(Y),
    FUN = constructor,
    Ydat = Y,
    Xdat = X,
    radius = radius,
    dimension = dimension,
    mc.cores = nc
  ),
  silent = T)

  if (is(ANS, 'try-error')) {
    #Windows does not support mclapply...
    HOMOLOGY <- lapply(
      X = 1:ncol(Y),
      FUN = constructor,
      Ydat = Y,
      Xdat = X,
      radius = radius,
      dimension = dimension,
    )
  } else {
    HOMOLOGY <<- ANS
  }

  return(HOMOLOGY)
} #end-optNP

normalizeX <- function (X,
                       MAXS = NULL,
                       MINS = NULL,
                       inv = F) {
  #DOC
  #Normalizes a vector or a matrix according to the given 'MAXS' and 'MINS'.
  #If 'MAXS' and 'MINS' are not provided 'X' is scaled so that each column is within [0,1].
  #ARGUMENTS
  #X: matrix or vector to normalize
  #MAXS, MINS: maxima and minima. NULL or NA values are set to max(X[,i]) and min(X[,i]) respectively.
  #inv: if T performs the inverse transformation
  #MAIN
  X = as.matrix(X)
  if (inv) {
    return(scale(
      scale(X, center = F, scale = (MAXS - MINS) ** (-1)),
      center = -MINS,
      scale = F
    ))
  } else {
    if (is.null(MAXS))
      MAXS = apply(X, 2, max)
    if (is.null(MINS))
      MINS = apply(X, 2, min)
    for (i in 1:ncol(X)) {
      if (is.na(MAXS[i]))
        MAXS[i] = max(X[, i])
      if (is.na(MINS[i]))
        MINS[i] = min(X[, i])
    }
    return(scale(X, center = MINS, scale = MAXS - MINS))
  }
}

plot.sobolManifold <- function(x, ...) {
  yrng = range(c(1 + x[['S']][, 'se'], x[['S']][, 'q0.05']))
  #plot estimates
  plot(
    x = 1:nrow(x[['S']]),
    y = x[['S']][, 'Si'],
    pch = 19,
    ylim = yrng,
    xaxt = 'n',
    xlab = 'parameter',
    ylab = 'Si and se',
    ...
  )
  #plot q 0.05
  points(
    x = 1:nrow(x[['S']]),
    y = x[['S']][, 'q0.05'],
    pch = 19,
    col = 2
  )
  #plot se
  arrows(
    x0 = 1:nrow(x[['S']]),
    y0 = x[['S']][, 'Si'],
    y1 = x[['S']][, 'Si'] - x[['S']][, 'se'],
    length = 0
  )	#lower
  arrows(
    x0 = 1:nrow(x[['S']]),
    y0 = x[['S']][, 'Si'],
    y1 = x[['S']][, 'Si'] + x[['S']][, 'se'],
    length = 0
  )	#upper
  #plot 0
  abline(h = 0, lty = 2)
  #x axis
  axis(side = 1,
       at = 1:nrow(x[['S']]),
       labels = row.names(x[['S']]))
  #legend
  legend(
    'topright',
    legend = c('Si', 'se', 'q0.05'),
    pch = c(19, NA, 19),
    lty = c(NA, 1, NA),
    col = c(1, 1, 2),
    horiz = T,
    bty = 'n'
  )
}

print.sobolManifold <- function(x, ...) {
  cat("\nCall:\n", deparse(x[['call']]), "\n", sep = "")
  cat("\nNumber of observations:", length(x[['Y']]), "\n")
  cat("\nFirst order indices:\n")
  print(x[['S']])
}

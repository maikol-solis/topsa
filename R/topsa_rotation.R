#topsa_angulos
#OJO: También hay que correr la función estimate_symmetric_reflection de helpers
topsa_angles <-
  function(Ydat,
           Xdat,
           #threshold.area = 0.9,
           threshold.radius = rep(0.05,ncol(Xdat)),
           # knearest = 20,
           method = "Alpha",angle) {
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
    # }#end-require-sf
    # # if (!requireNamespace("rgeos", quietly = TRUE)) {
    # #   stop("Please install the package rgeos: install.packages('rgeos')")
    # # }#end-require-rgeos
    # if (!requireNamespace("pbmcapply", quietly = TRUE)) {
    #   stop("Please install the package pbmcapply: install.packages('pbmcapply')")
    # }#end-require-pbmcapply
    # # if (!requireNamespace("multimode", quietly = TRUE)) {
    # #   stop("Please install the package multimode: install.packages('multimode')")
    # # }#end-require-multimode

    # future::plan("multisession")
    Xdat <- as.data.frame(Xdat)
    Ydat <- as.data.frame(Ydat)

    ANS <- list()
    ANS[['call']] <- match.call()
    ANS[['Xdat']] <- Xdat
    ANS[['Ydat']] <- Ydat
    # Xs <- lapply(X=1:NCOL(Xdat),function(i){scales::rescale(Xdat[,i])})
    # Ys <- scales::rescale(Ydat)
    # #Se centra para el calculo
    # meanXs <- unlist(lapply(X=1:NCOL(Xdat),function(i){mean(Xs[[i]])}))
    # meanYs <- mean(Ys)
    # XsYs <- lapply(X=1:ncol(Xdat),function(i){data.frame(Xs[[i]]-meanXs[i],Ys-meanYs)})
    # Rotated <- lapply(X=1:ncol(Xdat),function(i){as.matrix(XsYs[[i]]) %*% RotMat(-angles)})
    # #Se descentra para el calculo
    # Xs <- lapply(X=1:ncol(Xdat), function(i){Rotated[[i]][,1] + meanXs[i] })
    # Ys <- lapply(X=1:ncol(Xdat), function(i){Rotated[[i]][,2] + meanYs[i] })
    # ANS[['Xs']] <- dplyr::bind_cols(Xs)
    # ANS[['Ys']] <- Ys
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

    # cores <- parallel::detectCores(logical = FALSE)
    # cl <- parallel::makeCluster(cores)
    # parallel::clusterExport(cl,
    #                         c('topsa:::estimate_sensitivity_index',
    #                           'topsa:::VR_homology',
    #                           'topsa:::estimate_symmetric_reflection',
    #                           'Ydat',
    #                           'Xdat',
    #                           'threshold.radius',
    #                           'method'
    #                         ), envir=environment())


    sensitivity_results <- try(parallel::mclapply(
      X = 1:ncol(Xdat),
      FUN = estimate_sensitivity_index_angle,
      Ydat = Ydat,
      Xdat = Xdat,
      #   dimension = dimension,
      #   knearest = knearest,
      threshold = threshold.radius,
      method = method,
      angle=angle, #Se agrega angle
      mc.cores = parallel::detectCores(logical = FALSE)
    ),
    silent = T)



    if (is(sensitivity_results, 'try-error')){
      sensitivity_results <-
        lapply(
          # cl = cl,
          X = 1:ncol(Xdat),
          FUN = estimate_sensitivity_index_angle,
          Ydat = Ydat,
          Xdat = Xdat,
          # dimension = dimension,
          # knearest = knearest,
          # threshold.area = threshold.area,
          threshold = threshold.radius,
          method = method,
          angle=angle #Se agrega angle
        )
    }


    for (i in 1:ncol(Xdat)) {
      sensitivity_results[[i]]$xname <- colnames(Xdat)[i]
      sensitivity_results[[i]]$yname <- colnames(Ydat)
    }

    ANS[["results"]] <- sensitivity_results
    class(ANS) <- 'topsa'
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



estimate_sensitivity_index_angle <- function(ivar,
                                       Ydat,
                                       Xdat,
                                       dimension,
                                       # knearest,
                                       threshold, method,angle) {#Se agrega el parametro angle
  constructHOMOLOGY_angle <-
    function (ivar, Ydat, Xdat, dimension, threshold,method,angle) {#Se agrega el parametro angle
      Y <- as.matrix(Ydat)
      X <- as.matrix(Xdat[, ivar])

      idx <- order(X,Y)
      X <- X[idx,]
      Y<- Y[idx,]

      RotMat <- function(angle){
        matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), nrow=2, ncol=2)
      }

      #Aqui se hace el cambio!
      X <- as.matrix(scales::rescale(X , to = c(0, 1)), ncol = 1)
      Y <- as.matrix(scales::rescale(Y , to = c(0, 1)), ncol = 1)
      meanX <- mean(X)
      meanY <- mean(Y)
      #Se centra para el calculo
      XsYs <- data.frame(X-meanX,Y-meanY)
      Rotated <- as.matrix(XsYs) %*% RotMat(-angle)
      #Se descentra para el calculo
      X <- Rotated[,1] + meanX
      Y <- Rotated[,2] + meanY
      if (method == "Alpha") {
        Filtration <- TDA::alphaComplexFiltration(cbind(Y, X),printProgress = FALSE)
        cmplx <- Filtration$cmplx[Filtration$values<=threshold[ivar]]
      } else if(method=="VR"){
        Filtration <-TDA::ripsFiltration(
          cbind(Y, X),
          maxdimension = 1,
          maxscale = threshold[ivar],
          printProgress = FALSE
        )
        cmplx <- Filtration$cmplx
      }else{
        Filtration <- NULL
        stop("No method defined")
      }


      idx_triangles <- lengths(Filtration$cmplx) == 3
      clq <- cmplx[idx_triangles]
      # clq <- igraph::cliques(graphBase, min = 3, max = 3)
      clq <- matrix(unlist(clq), ncol = 3, byrow = TRUE)
      clq <- cbind(clq, clq[, 1])
      clq <- clq[order(clq[, 1], -clq[, 2], clq[, 3]), ]

      clq_list <-
        try(parallel::mclapply(
          X = seq_len(nrow(clq)),
          FUN =  function(i) {
            p <- sf::st_polygon(list(cbind(X[clq[i, ]], Y[clq[i, ]])))
          },
          mc.cores = parallel::detectCores(logical = FALSE)
        ),silent = T)


      if (is(clq_list, 'try-error')) {
        clq_list <-
          lapply(
            X = seq_len(nrow(clq)),
            FUN =  function(i) {
              p <- sf::st_polygon(list(cbind(X[clq[i, ]], Y[clq[i, ]])))
            }
          )
      }


      # message("Assemblying all the polygons and setting the manifold...")
      clq_polygons <-sf::st_geometrycollection(clq_list)

      clq_polygons <- sf::st_collection_extract(clq_polygons)

      grid_polygons <-
        sf::st_make_grid(clq_polygons, cellsize = 0.01)

      Number_Polygons <-
        lengths(sf::st_intersects(grid_polygons, clq_polygons))

      grid_polygons <-
        sf::st_sf(grid_polygons, Number_Polygons = Number_Polygons)

      clq_polygons <- sf::st_union(clq_polygons)

      #
      #         for (d in 1:dimension) {
      #           clq <- igraph::cliques(graphBase, min = d, max = d)
      #           clq <- matrix(unlist(clq),ncol = d,byrow = TRUE)
      #           if (d == 1) {
      #             H[[d]] <- st_multipoint(cbind(X[clq],Y[clq]))
      #             })
      #           } else if (d == 2) {
      #
      #             clq_list  <- purrr::map(clq, function(x) {
      #               cbind(X[x], Y[x])
      #             })
      #             H[[d]] <-  sf::st_multilinestring(clq_list)
      #           }
      #           else{
      #             clq_list  <- purrr::map(clq, function(x) {
      #               list(cbind(X[c(x, x[1])], Y[c(x, x[1])]))
      #             })
      #
      #             H[[d]] <-  sf::st_multipolygon(clq_list)
      #           }
      #         }





      return(
        list(
          grid_polygons = grid_polygons,
          manifold_unioned = clq_polygons,
          # neigborhood.distance = neigborhood.distance,
          threshold = threshold[ivar]
          # Number.Edges.per.Point = npositives
        )
      )


    } #end-function-constructor

  Ydat <- as.matrix(Ydat)

  H <- constructHOMOLOGY_angle(ivar, Ydat, Xdat, dimension, threshold, method,angle)

  # vv <- igraph::as_data_frame(H[["graph"]], "vertices")
  # H2 <- H[["homology"]][[3]]
  # l <- list()
  # idxObj <- sort(unique(as.numeric(H2)))
  # Vertices <- data.frame (name = vv$name[idxObj],
  #                         x = vv$x[idxObj],
  #                         y = vv$y[idxObj])
  #
  # l <- furrr::future_map(1:nrow(H2),
  #                        function(i) {
  #                          idxTriangle <- H2[i,]
  #                          idxname <-
  #                            Vertices$name %in% c(idxTriangle, idxTriangle[1])
  #                          subVertices <-  Vertices[idxname , ]
  #                          subVertices <-
  #                            rbind(subVertices, subVertices[1,])
  #                          Triangle <-
  #                            as.matrix((subVertices[, c(2, 3)]))
  #                          rownames(Triangle) <- NULL
  #                          p <- sf::st_polygon(list(Triangle))
  #                        }, .progress = TRUE)


  #     l <-  try(pbmcapply::pbmclapply(
  #   X = 1:nrow(H2),
  #   FUN = function(i) {
  #     idxTriangle <- H2[i, ]
  #     idxname <- Vertices$name %in% c(idxTriangle, idxTriangle[1])
  #     subVertices <-  Vertices[idxname ,]
  #     subVertices <- rbind(subVertices, subVertices[1, ])
  #     Triangle <- as.matrix((subVertices[, c(2, 3)]))
  #     rownames(Triangle) <- NULL
  #     p <- sf::st_polygon(list(Triangle))
  #   },
  #   mc.cores = mc.cores.inner
  # ),
  # silent = T)

  # if (is(l, 'try-error')) {
  #   #Windows does not support mclapply...
  #   l <- lapply(
  #     X = 1:nrow(H2),
  #     FUN = function(i) {
  #       idxTriangle <- H2[i, ]
  #       idxname <- Vertices$name %in% c(idxTriangle, idxTriangle[1])
  #       subVertices <-  Vertices[idxname ,]
  #       subVertices <- rbind(subVertices, subVertices[1, ])
  #       Triangle <- as.matrix((subVertices[, c(2, 3)]))
  #       rownames(Triangle) <- NULL
  #       p <- sf::st_polygon(list(Triangle))
  #     }
  #   )
  # }

  # mp <- sf::st_multipolygon(l)
  #mp_union <- sf::st_buffer(sf::st_union(mp),0)
  #

  # mp_union <- sf::st_simplify(H[["homology"]][[3]])

  mp_union <- H[["manifold_unioned"]]

  mp_reflection <- estimate_symmetric_reflection(mp_union)

  bb <- sf::st_make_grid(x = mp_union, n = 1)

  mp_sym_difference <-
    sf::st_sym_difference(mp_union, mp_reflection)


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


library(ggplot2)
#Función data.frame-----
topsa_index_angle <- function(Ydat,
                             Xdat,
                             #threshold.area = 0.9,
                             threshold.radius = rep(0.05,ncol(Xdat)),
                             # knearest = 20,
                             method = "Alpha",steps=pi/6){ #Por el momento Xdat solo será una columna



  #cita <- seq(0,2*pi,steps) # Se elimina porque serepiten
  Xdat <- as.matrix(Xdat)
  Ydat <- as.matrix(Ydat)
  cita <- seq(0,pi,steps)
  topsa_by_angle_col <- list()

  for (j in 1:ncol(Xdat)) {
    topsa_by_angle_col[[j]] <- data.frame(Xcol=j,Index=unlist(lapply(X=1:length(cita), function(i){estimate_sensitivity_index_angle(j,Ydat,
                                                                           Xdat,
                                                                           dimension,
                                                                           threshold = threshold.radius,
                                                                           method,angle=cita[i])$Symmetric.Diff.Area})))
  }
  data_topsa_by_angle_col <- dplyr::bind_rows(topsa_by_angle_col)
  data_topsa_by_angle_col[,"Angle"] <- cita
  data_topsa_by_angle_col
}
#Función plot-----
topsa_plot_angle <- function(Ydat,
                             Xdat,
                             #threshold.area = 0.9,
                             threshold.radius = rep(0.05,ncol(Xdat)),
                             # knearest = 20,
                             method = "Alpha",steps=pi/6){ #Por el momento Xdat solo será una columna



  #cita <- seq(0,2*pi,steps) # Se elimina porque serepiten
  Xdat <- as.matrix(Xdat)
  Ydat <- as.matrix(Ydat)
  cita <- seq(0,pi,steps)
  topsa_by_angle_col <- list()

  for (j in 1:ncol(Xdat)) {
    topsa_by_angle_col[[j]] <- data.frame(Xcol=j,Index=unlist(lapply(X=1:length(cita), function(i){estimate_sensitivity_index_angle(j,Ydat,
                                                                                                                                    Xdat,
                                                                                                                                    dimension,
                                                                                                                                    threshold = threshold.radius,
                                                                                                                                    method,angle=cita[i])$Symmetric.Diff.Area})))
  }
  data_topsa_by_angle_col <- dplyr::bind_rows(topsa_by_angle_col)
  data_topsa_by_angle_col[,"Angle"] <- cita
  p <- ggplot(data_topsa_by_angle_col) + geom_line(aes(Angle,Index)) + geom_point(aes(Angle,Index)) + facet_wrap(~Xcol)
  # Agregar como helper
  fracAx <- function(p, symbol="pi", width=steps/pi) {
    require(MASS)                                                 # for fractions
    val <- tryCatch(eval(parse(text=symbol)), error=function(e) 1)
    info <- ggplot_build(p)
    xrange <- c(0,pi)                 # get the x-range of figure
    vec.breaks <- seq(floor(xrange[1]), ceiling(xrange[2]), by=width)
    fracs <- strsplit(attr(fractions(vec.breaks), "fracs"), "/")  # convert to fractions
    labels <- sapply(fracs, function(i)
      if (length(i) > 1) { paste(i[1], "*", symbol, "/", i[2]) }
      else { paste(i, "*", symbol) })
    p + scale_x_continuous(breaks=vec.breaks*val, labels=parse(text=labels))
  }
  ## Make the graph with pi axis
  fracAx(p, "pi")

}

#Pruebas -----












ishigami.fun <- function(X) {
  A <- 7
  B <- 0.1
  sin(X[, 1]) + A * sin(X[, 2])^2 + B * X[, 3]^4 * sin(X[, 1])
}

X <- matrix(runif(3*1000, -pi, pi), ncol = 3)
Y <- ishigami.fun(X)


(p_plot <- topsa_plot_angle(Y,X,steps = pi/12,threshold.radius = c(0.003,0.002,0.001)))
(p_data <- topsa_index_angle(Y,X,steps = pi/12,threshold.radius = c(0.003,0.002,0.001)))




asdh <- topsa_angles(Y,X,threshold.radius = (0.002),angle = 3*pi/4)
plot(asdh$results[[1]]$homology$manifold_unioned,axes=T)
topsa::print(asdh)
# p <- topsa_plot_angle(y,x,steps = pi/6)
plot(asd$results[[1]]$homology$manifold_unioned,axes=T)


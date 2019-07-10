#' plot \code{topsa} objects
#'
#' Plot method for objects of class \code{topsa}.
#'
#' @param topsaObj an object of class \code{topsa}
#' @param nvar  it could be a secuence from 1 to the number of variables
#'   indicating which variables should be plotted. It could be the character
#'   'all' for plot all the variables.
#' @param ... further arguments passed to the \code{plot} function
#'
#' @return A plot of generated with the output of \code{topsa}. For each
#' variable in the model, it creates the plot of the correponding manifold, its
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
#' estimation <- topsa(Ydat = Y, Xdat = X)
#'
#' plot(estimation)


plot <- function(topsaObj, nvar = "all", ...) {
  UseMethod("plot", topsaObj)
}

#' @export
#' @rdname plot
plot.topsa <- function(topsaObj,
                       nvar = "all",
                       ...) {
  if (class(topsaObj)!= "topsa"){
    stop("This function only accepts 'topsa' class objects")
  }

  if (is.character(nvar)) {
    if (nvar == "all") {
      nvar = seq(1, length(topsaObj$results))
    } else  {
      message("Please enter 'all' or vector of values to plot.")
    }
  } else if (!is.numeric(nvar)) {
    message("Please enter 'all' or vector of values to plot.")
  }




  manifold.all <- NULL

  for (k in nvar) {
    manifold <- topsaObj$result[[k]]$manifold.plot

    mp_origin_coords <- sf::st_coordinates(manifold)

    manifold <-
      manifold - c(min(mp_origin_coords[,"X"]), min(mp_origin_coords[,"Y"]))
    mp_coordinates <- sf::st_coordinates(manifold)

    bb <- sf::st_make_grid(x = manifold, n = 1)

    reflection <-  matrix(c(1, 0, 0, -1), 2, 2)

    #    xy_original <- extract_XY_coords(mp_union)

    manifold_reflectionx <- (manifold * reflection)

    # mp_coordinates_reflection <- sf::st_coordinates(mp_reflection)

    manifold_reflectionx <- manifold_reflectionx +
      c(0, 2 * mean(mp_coordinates[, "Y"]))
    # c(0, min(mp_coordinates_reflection[, "Y"]))


    boundingbox <- sf::st_make_grid(x = manifold, n = 1)

    datapoints <-
      sf::st_multipoint(as.matrix(cbind(topsaObj$Xdat[, k], topsaObj$Ydat)))





#
#     reflectionx <-  matrix(c(1, 0, 0, -1), 2, 2)
#
#     manifold_reflectionx <-
#       (manifold - c(min(xy_original$x), min(xy_original$y))) * reflectionx
#
#     xy_reflection <- extract_XY_coords(manifold_reflectionx)
#
#     manifold_reflectionx <-
#       manifold_reflectionx - c(0, min(xy_reflection$y) + max(xy_reflection$y))
#
#
#     manifold_reflectionx <-
#       manifold_reflectionx + c(min(xy_original$x), min(xy_original$y))
#
#
#     manifold_sym_difference <-
#       sf::st_sym_difference(manifold, manifold_reflectionx)

    manifold_intersection <-
      sf::st_intersection(manifold, manifold_reflectionx)

    manifold <-
      manifold + c(min(mp_origin_coords[, "X"]), min(mp_origin_coords[, "Y"]))
    manifold_reflectionx <-
      manifold_reflectionx + c(min(mp_origin_coords[, "X"]), min(mp_origin_coords[, "Y"]))
    manifold_intersection <-
      manifold_intersection + c(min(mp_origin_coords[, "X"]), min(mp_origin_coords[, "Y"]))
boundingbox <- boundingbox + c(min(mp_origin_coords[, "X"]), min(mp_origin_coords[, "Y"]))

    object_type <- NULL

    manifold.all <-
      rbind(
        manifold.all,
        sf::st_sf(
          variable = k,
          alpha = 1,
          object_type = "Bounding Box",
          geom = sf::st_geometry(boundingbox)
        ),
        sf::st_sf(
          variable = k,
          alpha = 1,
          object_type = "Original",
          geom = sf::st_geometry(manifold)
        ),
        sf::st_sf(
          variable = k,
          alpha = 1,
          object_type = "Sym. Reflection",
          geom = sf::st_geometry(manifold_reflectionx)
        ),
        sf::st_sf(
          variable = k,
          alpha = 0,
          object_type = "Intersection",
          geom = sf::st_geometry(manifold_intersection)
        ),
        sf::st_sf(
          variable = k,
          alpha = 1,
          object_type = "DataPoints",
          geom = sf::st_geometry(datapoints)
        )
      )
  }

  ggplot2::ggplot() +
    ggplot2::geom_sf(data = manifold.all,
                     ggplot2::aes(
                       fill = object_type,
                       size = object_type,
                       color = object_type
                     )) +
    ggplot2::scale_fill_manual(
      name = "Estimated Areas",
      values = c(
        "Original" = "#a6cee3",
        "Sym. Reflection" =  "#b2df8a",
        "Intersection" = "white",
        "Bounding Box" = NA,
        "DataPoints" = NA
      ),
      breaks = c("Original" ,
                 "Sym. Reflection",
                 "Intersection",
                 NA,
                 NA)
    ) +
    ggplot2::scale_size_manual(
      name = "",
      values = c(
        "Original" = 1,
        "Sym. Reflection" =  1,
        "Intersection" = 1,
        "Bounding Box" = 2,
        "DataPoints" = 0.1
      ),
      guide = FALSE
    ) +
    ggplot2::scale_color_manual(
      name = "",
      values = c(
        "Original" = "#1f78b4",
        "Sym. Reflection" = "#33a02c",
        "Intersection" = NA,
        "Bounding Box" = "#e31a1c",
        "DataPoints" = "grey50"
      ),
      guide = FALSE
    ) +
    ggplot2::facet_wrap( ~ variable) +
    ggplot2::coord_sf() +
    ggplot2::theme_minimal(...)
}

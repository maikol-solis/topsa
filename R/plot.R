#' Plot method for objects \code{TopSA}
#'
#' Plot the Complex, symmetric reflection and symmetric diffrence for and object
#' of class \code{\ĺink{TopSA}}
#'
#' @param TopSAObj an object of class \code{TopSA}
#' @param nvar  it could be a secuence from 1 to the number of variables
#'   indicating which variables should be plotted. It could be the character
#'   'all' for plot all the variables.
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

plot <- function(TopSAObj, nvar = "all", ...) {
  UseMethod("plot", TopSAObj)
}

# #' @export
# #' @rdname plot
plot.TopSA <- function(TopSAObj,
                       nvar = "all",
                       ...) {

  if (is.character(nvar)) {
    if (nvar == "all") {
      nvar = seq(1, length(TopSAObj$results))
    } else  {
      message("Please enter 'all' or vector of values to plot.")
    }
  } else if (!is.numeric(nvar)) {
    message("Please enter 'all' or vector of values to plot.")
  }




  manifold.all <- NULL

  for (k in nvar){



    manifold <- TopSAObj$result[[k]]$manifold.plot

    Xcoord <- manifold[[1]][, 1]
    Ycoord <- manifold[[1]][, 2]


    boundingbox <- sf::st_make_grid(x = manifold, n = 1)

    datapoints <-
      sf::st_multipoint(as.matrix(cbind(TopSAObj$Xdat[, k], TopSAObj$Ydat)))


    reflectionx <-  matrix(c(1, 0, 0, -1), 2, 2)

    manifold_reflectionx <-
      ((manifold - c(min(Xcoord), min(Ycoord))) * reflectionx +
         c(0,   min(Ycoord) +  max(Ycoord))) +
      c(min(Xcoord), min(Ycoord))


    manifold_sym_difference <-
      sf::st_sym_difference(manifold, manifold_reflectionx)

    manifold_intersection <-
      sf::st_intersection(manifold, manifold_reflectionx)

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
    ggplot2::geom_sf(data = manifold.all, ggplot2::aes(
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
      ), guide = FALSE
    ) +
    ggplot2::scale_color_manual(
      name = "",
      values = c(
        "Original" = "#1f78b4",
        "Sym. Reflection" = "#33a02c",
        "Intersection" = NA,
        "Bounding Box" = "#e31a1c",
        "DataPoints" = "grey50"
      ), guide = FALSE
    ) +
    ggplot2::facet_wrap(~variable) +
    ggplot2::coord_sf() +
    ggplot2::theme_minimal()
}


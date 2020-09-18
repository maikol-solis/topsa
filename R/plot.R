#' plot \code{topsa} objects
#'
#' Plot method for objects of class \code{topsa}.
#'
#' @param topsaObj an object of class \code{topsa}
#' @param nvar  it could be a sequence from 1 to the number of variables
#'   indicating which variables should be plotted. It could be the character
#'   'all' for plot all the variables.
#' @param ... further arguments passed to the \code{plot} function
#'
#' @return A plot of generated with the output of \code{topsa}. For each
#' variable in the model, it creates the plot of the corresponding manifold, its
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
#' X <- matrix(runif(3*50, -pi, pi), ncol = 3)
#' Y <- ishigami.fun(X)
#'
#' estimation <- topsa(Ydat = Y, Xdat = X)
#'
#' plot_topsa(estimation)

plot_topsa <- function(topsaObj,
                          nvar = "all",
                          ...) {
  if (class(topsaObj) != "topsa") {
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
    manifold <- topsaObj$result[[k]]$homology$manifold_unioned

    boundingbox <- sf::st_make_grid(x = manifold, n = 1)

    datapoints <-
      sf::st_multipoint(as.matrix(cbind(topsaObj$Xr[, k], topsaObj$Yr)))

    object_type <- NULL

    manifold.all <-
      rbind(
        manifold.all,
        sf::st_sf(
          variable = colnames(topsaObj$Xdat)[k],
          alpha = 1,
          object_type = "Bounding Box",
          geom = sf::st_geometry(boundingbox)
        ),
        sf::st_sf(
          variable = colnames(topsaObj$Xdat)[k],
          alpha = 1,
          object_type = "Manifold",
          geom = sf::st_geometry(manifold)
        ),
        sf::st_sf(
          variable = colnames(topsaObj$Xdat)[k],
          alpha = 0.5,
          object_type = "DataPoints",
          geom = sf::st_geometry(datapoints)
        )
      )

  }


  plotManifold <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = manifold.all,
      ggplot2::aes(fill = object_type,
                   # size = object_type,
                   color = object_type),
      shape = 20
    ) +
    ggplot2::scale_fill_manual(
      name = "Estimated Areas",
      values = c(
        "Manifold" = "#a6cee3",
        "Sym. Reflection" =  "#1f78b4",
        "Intersection" = "#b2df8a",
        "Bounding Box" = NA,
        "DataPoints" = NA
      ),
      breaks = c("Manifold" ,
                 "Sym. Reflection",
                 "Intersection",
                 NA,
                 NA)
    ) +
    ggplot2::scale_size_manual(
      name = "",
      values = c(
        "Manifold" = 1,
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
        "Manifold" = "#1f78b4",
        "Sym. Reflection" = "#a6cee3",
        "Intersection" = "#33a02c",
        "Bounding Box" = "#e31a1c",
        "DataPoints" = ggplot2::alpha("black",0.4)
      ),
      guide = FALSE
    ) +
    ggplot2::facet_wrap( ~ variable, ncol = 3) +
    ggplot2::coord_sf() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      aspect.ratio = 1,
      # panel.spacing = ggplot2::unit(1, "lines"),
      strip.background = ggplot2::element_rect(fill = "grey95"),
      # strip.placement = "outside"
    )


  return(plotManifold)
}

plot_sym_reflection <- function(topsaObj,
                                nvar = "all",
                                ...) {
  if (class(topsaObj) != "topsa") {
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




  manifold.sym.reflection.all <- NULL

  for (k in nvar) {
    manifold <- topsaObj$result[[k]]$homology$manifold_unioned

    manifold_reflectionx <- estimate_symmetric_reflection(manifold)

    boundingbox <- sf::st_make_grid(x = manifold, n = 1)

    datapoints <-
      sf::st_multipoint(as.matrix(cbind(topsaObj$Xr[, k], topsaObj$Yr)))

    manifold_intersection <-
      sf::st_intersection(manifold, manifold_reflectionx)

    object_type <- NULL

    manifold.sym.reflection.all <-
      rbind(
        manifold.sym.reflection.all,
        sf::st_sf(
          variable = colnames(topsaObj$Xdat)[k],
          alpha_var = 1,
          object_type = "Bounding Box",
          geom = sf::st_geometry(boundingbox)
        ),
        sf::st_sf(
          variable = colnames(topsaObj$Xdat)[k],
          alpha_var = 1,
          object_type = "Manifold",
          geom = sf::st_geometry(manifold)
        ),
        sf::st_sf(
          variable = colnames(topsaObj$Xdat)[k],
          alpha_var = 1,
          object_type = "Sym. Reflection",
          geom = sf::st_geometry(manifold_reflectionx)
        ),
        sf::st_sf(
          variable = colnames(topsaObj$Xdat)[k],
          alpha_var = 0,
          object_type = "Intersection",
          geom = sf::st_geometry(manifold_intersection)
        ),
        sf::st_sf(
          variable = colnames(topsaObj$Xdat)[k],
          alpha_var = 0.5,
          object_type = "DataPoints",
          geom = sf::st_geometry(datapoints)
        )
      )

  }


  plotSymmetricReflection <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = manifold.sym.reflection.all,
      ggplot2::aes(fill = object_type,
                   color = object_type),
      shape = 20
    ) +
    ggplot2::facet_wrap( ~ variable, ncol = 3) +
    ggplot2::scale_fill_manual(
      name = "Estimated Areas",
      values = c(
        "Manifold" = "#a6cee3",
        "Sym. Reflection" =  "#1f78b4",
        "Intersection" = "#b2df8a",
        "Bounding Box" = NA,
        "DataPoints" = NA
      ),
      breaks = c("Manifold" ,
                 "Sym. Reflection",
                 "Intersection",
                 NA,
                 NA)
    ) +
    ggplot2::scale_size_manual(
      name = "",
      values = c(
        "Bounding Box" = 2,
        "Manifold" = 1,
        "Sym. Reflection" =  1,
        "Intersection" = 1,
        "DataPoints" = NA
      ),
      guide = FALSE
    ) +
    ggplot2::scale_color_manual(
      name = "",
      values = c(
        "Manifold" = "#1f78b4",
        "Sym. Reflection" = "#a6cee3",
        "Intersection" = "#33a02c",
        "Bounding Box" = "#e31a1c",
        "DataPoints" = ggplot2::alpha("black",0.4)
      ),
      guide = FALSE
    ) +
    # ggplot2::facet_wrap( ~ variable, ncol = 3) +
    ggplot2::coord_sf()  +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      aspect.ratio = 1,
      # panel.spacing = ggplot2::unit(1, "lines"),
      strip.background = ggplot2::element_rect(fill = "grey95"),
    )


  # plotIntersection <-  ggplot2::ggplot() +
  #   ggplot2::geom_sf(data = manifold.sym.difference.all,
  #                    mapping = ggplot2::aes(
  #                      fill = object_type,
  #                      size = object_type,
  #                      color = object_type
  #                    )) +
  #   ggplot2::scale_fill_manual(
  #     name = "Estimated Areas",
  #     values = c(
  #       "Original" = "#a6cee3",
  #       "Sym. Reflection" =  "#b2df8a",
  #       "Intersection" = "orange",
  #       "Bounding Box" = NA,
  #       "DataPoints" = NA
  #     ),
  #     breaks = c("Original" ,
  #                "Sym. Reflection",
  #                "Intersection",
  #                NA,
  #                NA)
  #   ) +
  #   ggplot2::scale_size_manual(
  #     name = "",
  #     values = c(
  #       "Original" = 1,
  #       "Sym. Reflection" =  1,
  #       "Intersection" = 1,
  #       "Bounding Box" = 2,
  #       "DataPoints" = 0.1
  #     ),
  #     guide = FALSE
  #   ) +
  #   ggplot2::scale_color_manual(
  #     name = "",
  #     values = c(
  #       "Original" = "#1f78b4",
  #       "Sym. Reflection" = "#33a02c",
  #       "Intersection" = "darkorange",
  #       "Bounding Box" = "#e31a1c",
  #       "DataPoints" = "grey50"
  #     ),
  #     guide = FALSE
  #   ) +
  #   ggplot2::facet_wrap( ~ variable) +
  #   ggplot2::coord_sf()  +
  #   ggplot2::theme_minimal() +
  #   ggplot2::theme(panel.spacing = ggplot2::unit(1, "lines"))

  return(plotSymmetricReflection)
}

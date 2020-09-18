#' Barcode plotter
#'
#' @param Ydat A vector with the model's dependent variable.
#' @param Xdat A matrix with the model's input variables.
#' @param maxscale Maximum radius allowed to find the barcode.
#' @param mc.cores Number of cores used to estimate the barcodes in parallel. (See \code{\link{mclapply}}).
#'
#' @return A plot with the barcode for each variable.
#' @export
#'
#' @examples
#'ishigami.fun <- function(X) {
#' A <- 7
#' B <- 0.1
#' sin(X[, 1]) + A * sin(X[, 2])^2 + B * X[, 3]^4 * sin(X[, 1])
#' }

#' X <- matrix(runif(3*100, -pi, pi), ncol = 3)
#' Y <- ishigami.fun(X)
#'
#'barcode_plotter(Ydat = Y, Xdat = X, maxscale = 0.2, mc.cores = 2)
#'
#' @importFrom stats median quantile
#' @importFrom utils head
#' @importFrom methods is

barcode_plotter <- function(Ydat,
                            Xdat,
                            maxscale = rep(0.05, ncol(Xdat)),
                            mc.cores = parallel::detectCores(logical = FALSE)) {
  if (length(maxscale) == 1) {
    maxscale <- rep(maxscale, ncol(Xdat))
  } else if (length(maxscale) < ncol(Xdat)) {
    stop("A vector of size 1 or ncol(Xdat) with numbers between 0 or 1 must be entered.")
  }

  Xdat <- as.data.frame(Xdat)
  Ydat <- as.data.frame(Ydat)

  Xr <- as.data.frame(lapply(Xdat, scales::rescale))
  Yr <- as.data.frame(lapply(Ydat, scales::rescale))




  p <-
    try(parallel::mclapply(
      X = seq_len(ncol(Xdat)),
      FUN = function(i) {
        d <- TDA::ripsDiag(
          cbind(Yr, Xr[, i]),
          maxdimension = 1,
          maxscale = maxscale[i],
          library = "GUDHI",
          printProgress = FALSE
        )


        d$diagram <- d$diagram[-1,]

        ord <-
          order(d$diagram[, "Death"] - d$diagram[, "Birth"], decreasing = TRUE)

        dim1 <- d$diagram[ord, "dimension"] == 1
        d$diagram[head(ord[dim1]), "dimension"] <- -1

        dim1 <- d$diagram[, "dimension"] == 1
        dim0 <- d$diagram[, "dimension"] == 0

        list(
          df_diagram = data.frame(
            variable = names(Xr)[i],
            d$diagram,
            height = 1:nrow(d$diagram)
          ),
          df_labels = data.frame(
            variable = names(Xr)[i],
            max_death0 = max(d$diagram[dim0, "Death"]),
            max_death1 = max(d$diagram[dim1 |
                                         d$diagram[, "dimension"] == -1, "Death"]),
            median_death1 = median(d$diagram[dim1, "Death"]),
            mean_height = mean(1:nrow(d$diagram))
          )
        )
      },
      mc.cores = mc.cores
    ),silent = T)

  if(is(p,'try-error')){
    p <- lapply(X = seq_len(ncol(Xdat)),
      FUN = function(i) {
        d <- TDA::ripsDiag(
          cbind(Yr, Xr[, i]),
          maxdimension = 1,
          maxscale = maxscale[i],
          library = "GUDHI",
          printProgress = FALSE
        )


        d$diagram <- d$diagram[-1,]

        ord <-
          order(d$diagram[, "Death"] - d$diagram[, "Birth"], decreasing = TRUE)

        dim1 <- d$diagram[ord, "dimension"] == 1
        d$diagram[head(ord[dim1]), "dimension"] <- -1

        dim1 <- d$diagram[, "dimension"] == 1
        dim0 <- d$diagram[, "dimension"] == 0

        list(
          df_diagram = data.frame(
            variable = names(Xr)[i],
            d$diagram,
            height = 1:nrow(d$diagram)
          ),
          df_labels = data.frame(
            variable = names(Xr)[i],
            max_death0 = max(d$diagram[dim0, "Death"]),
            max_death1 = max(d$diagram[dim1 |
                                         d$diagram[, "dimension"] == -1, "Death"]),
            median_death1 = median(d$diagram[dim1, "Death"]),
            mean_height = mean(1:nrow(d$diagram))
          )
        )
      }
    )
  }




  df_diagram <- lapply(p, function(x)
    x$df_diagram)

  df_labels <- lapply(p, function(x)
    x$df_labels)

  df_diagram <- do.call(rbind, df_diagram)
  df_labels <- do.call(rbind, df_labels)


  print(df_labels)


  ggplot2::ggplot(df_diagram) +
    ggplot2::geom_segment(
      mapping = ggplot2::aes(
        x = Birth,
        xend = Death,
        y = height,
        yend = height,
        color = as.factor(dimension)
      ),
      size = 1
    ) +
    ggplot2::facet_wrap(~ variable, scales = "free") +
    ggplot2::scale_color_manual(
      "Features",
      values = c("0" = "black", "1" = "dodgerblue3", "-1" = "tomato1"),
      limits = c("-1", "1", "0"),
      labels = c("Longest H1", "H1", "H0")
    ) +
    ggplot2::geom_vline(
      data = df_labels,
      ggplot2::aes(xintercept = median_death1),
      col = "dodgerblue3",
      linetype = "dashed"
    ) +
    ggplot2::geom_vline(
      data = df_labels,
      ggplot2::aes(xintercept = max_death1),
      col = "tomato1",
      linetype = "dashed"
    ) +
    ggplot2::geom_vline(
      data = df_labels,
      ggplot2::aes(xintercept = max_death0),
      col = "black",
      linetype = "dashed"
    ) +
    ggplot2::geom_label(
      data = df_labels,
      ggplot2::aes(
        label = paste0("H1: ", round(median_death1, 3)),
        x = median_death1,
        y = quantile(df_diagram$height, 0.5),
      ),
      color = "dodgerblue4",
      fill = "grey95",
      vjust = "inward",
      hjust = "inward"
    ) +

    ggplot2::geom_label(
      data = df_labels,
      ggplot2::aes(
        label = paste0("H1: ", round(max_death1, 3)),
        x = max_death1,
        y = max(df_diagram$height),
      ),
      color = "tomato3",
      fill = "grey95",
      vjust = "inward",
      hjust = "inward"
    ) +

    ggplot2::geom_label(
      data = df_labels,
      ggplot2::aes(
        label = paste0("H0: ", round(max_death0, 3)),
        x = max_death0,
        y = quantile(df_diagram$height, 0.1),
      ),
      fill = "grey95",
      vjust = "inward",
      hjust = "inward"
    ) +
    # ggplot2::geom_vline(data = df_labels,
    #                     ggplot2::aes(xintercept = q75_death),
    #                     col = "darkgreen") +
    # ggplot2::geom_text(
    #   data = df_labels,
    #   ggplot2::aes(
    #     label = "Q 75%",
    #     x = q75_death ,
    #     y = mean_height,
    #     angle = 90
    #   ),
  #   vjust = -0.5
  # ) +
  ggplot2::theme_classic() +
    ggplot2::theme(
      aspect.ratio = 1,
      #legend.position = "none",
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::ylab("") +
    ggplot2::xlab("Radius")

  # invisible(df_diagram)
}

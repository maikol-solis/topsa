
topsa_full_rotation <-  function(Ydat,
                                 Xdat,
                                 threshold.radius = rep(0.05, ncol(Xdat)),
                                 method = "Alpha",
                                 step =pi/6) {
  angles <- seq(0, pi-step, step)
  topsa_list_angles <-
    lapply(
      X = angles,
      FUN = function(a) {
        topsa(Ydat,
              Xdat,
              threshold.radius,
              method = "Alpha",
              angle = a)
      }
    )

  geometric_list_angles <- lapply(topsa_list_angles, function(x) {
    data.frame(
      variable = colnames(x$Xdat),
      index = sapply(x$results, function(y)
        y$Geometric.Sensitivity),
      angle = x$angle
    )
  })

  geometric_df_angles <-
    do.call(what = "rbind", args = geometric_list_angles)

  return(geometric_df_angles)
}

topsa_full_rotation_plot <- function(topsaRotationDF) {
  p <- ggplot(topsaRotationDF) +
    geom_line(aes(angle, index)) +
    geom_point(aes(angle, index)) +
    facet_wrap( ~ variable)

  return(p)
}


topsa_full_rotation_ecdf <- function(topsaRotationDF) {
  # topsaRotationDF_by <-    by(
  #   data = topsaRotationDF,
  #   INDICES = topsaRotationDF$variable,
  #   FUN = function(x) {
  #     data.frame(
  #       variable = unique(x$variable),
  #       min_index  = min(x$index),
  #       max_index = max(x$index)
  #     )
  #   }
  # )
  #
  # topsaRotationDF_by <- do.call("rbind", topsaRotationDF_by)
  #
  # topsaRotationDF <- merge(topsaRotationDF, topsaRotationDF_by, by = "variable")

  p <- ggplot(topsaRotationDF) +
    stat_ecdf(aes(index)) +
    stat_function(mapping = aes(index),
                  fun = punif,
                  # args =list(min = min_index, max = max_index),
                  color = "red") +
    facet_wrap(~ variable)
  return(p)
}

topsa_full_rotation_stats <- function(topsaRotationDF) {
  topsaRotationDF %>%
    group_by(variable) %>%
    summarise(sigma = sd(index),
              mu = mean(index),
              cv = sigma / mu)
}



topsa_full_rotation_test <- function(topsaRotationDF) {
  #usar otras pruebaa¿?
  variables <- as.vector(unique(topsaRotationDF$variable))
  results <- lapply(variables, function(i){
    vect <- topsaRotationDF %>% dplyr::filter(variable==i) %>% dplyr::select(index)
    ks.test(vect,"punif",min(vect),max(vect))
  })
  names(results) <- variables
  results
}
#topsa_full_rotation_test(estimation)









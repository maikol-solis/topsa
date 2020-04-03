estimate_symmetric_reflection <- function(polygon) {
  true_polygon_coords <- sf::st_coordinates(polygon)

  polygon <-
    polygon - c(min(true_polygon_coords[, "X"]), min(true_polygon_coords[, "Y"]))

  polygon_coords <- sf::st_coordinates(polygon)

  affine_transformation <-  matrix(c(1, 0, 0, -1), 2, 2)

  polygon_reflected <- (polygon * affine_transformation)
  polygon_reflected <- polygon_reflected +
    c(0, min(polygon_coords[, "Y"]) + max(polygon_coords[, "Y"]))
    #c(0, 2 * mean(polygon_coords[, "Y"]))
    #c(0, 2 * diff(range(polygon_coords[, "Y"])) / 2)

  polygon_reflected + c(min(true_polygon_coords[, "X"]), min(true_polygon_coords[, "Y"]))
}






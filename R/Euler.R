# Euler

# ishigami.fun <- function(X) {
#   A <- 7
#   B <- 0.1
#   sin(X[, 1]) + A * sin(X[, 2])^2 + B * X[, 3]^4 * sin(X[, 1])
# }
#
# X <- matrix(runif(3*100, -pi, pi), ncol = 3)
# Y <- ishigami.fun(X)
# Xdat <- X
# Ydat <- Y

#Funcion para seccionar los datos por angulo ----
#
Seccion_angulos <-
  function(Yr, Xr, steps = pi / 6) {
    # Reciben vector Y, X_i y el paso
    X_i <- Xr - mean(Xr)
    Y_i <- Yr - mean(Yr)
    angulos <- atan2(Y_i, X_i)


    angulos_data <- list(NULL)
    angulos <- ifelse(angulos >= 0, angulos, angulos + 2 * pi)
    cita <- seq(0, 2 * pi - steps, by = steps)
    for (l in 1:length(cita)) {
      if (cita[l] <= pi) {
        angulos_data[[l]] <-
          cbind(x = X_i[which(angulos > cita[l] &
                                angulos < cita[l] + pi)],
                y = Y_i[which(angulos > cita[l] &
                                angulos < cita[l] + pi)],
                angle = cita[l])
      } else{
        angulos_data[[l]] <-
          cbind(x = X_i[which(angulos >= cita[l] |
                                angulos <= -pi + cita[l])],
                y = Y_i[which(angulos >= cita[l] |
                                angulos <= -pi + cita[l])],
                angle = cita[l])

      }

    }
    names(angulos_data) <- cita

      # paste("cita", 0:(length(cita) - 1), sep = "")

    angulos_data
  }

pru <- Seccion_angulos(Y, X[, 1], pi / 2)
ggplot() + geom_point(aes(X[, 1] - mean(X[, 1]), Y - mean(Y)), shape = 2) + geom_point(aes(pru[[3]][, 1], pru[[3]][, 2]))


##### Característica de Euler ----
### Observación: Se solicitó con AlphaShapes, pero recordar que a esta función no se le puede establecer el maxscale

caracteristica <- function(Ydat,
                            Xdat,
                            maxscale = 0.05) {
  Xdat <- as.data.frame(Xdat)
  Ydat <- as.data.frame(Ydat)

  Xr <- as.data.frame(lapply(Xdat, scales::rescale))
  Yr <- as.data.frame(lapply(Ydat, scales::rescale))
  d <- TDA::alphaComplexDiag(cbind(Yr, Xr),
                         maxdimension = 1,
                         library = "GUDHI",
                         printProgress = FALSE)
  # d <- TDA::ripsDiag(
  #               cbind(Yr, Xr),
  #               maxdimension = 1,
  #               maxscale = maxscale,
  #               library = "GUDHI",
  #               printProgress = FALSE)


    d$diagram <- d$diagram[-1,]

    # ord <-
    #   order(d$diagram[, "Death"] - d$diagram[, "Birth"], decreasing = TRUE)
    #
    # dim1 <- d$diagram[ord, "dimension"] == 1
    # d$diagram[head(ord[dim1]), "dimension"] <- -1
d$diagram
}
#caracteristica(Y,X[,1],0.2)


#Esta función recibe la columna  Y, X_i
X_carac <- function(Y,X,maxscale) {
  Carac <- as.data.frame(caracteristica(Ydat = Y, Xdat = X, maxscale))
  #Carac <- Carac[Carac$dimension!=-1,]
  H0 <- Carac[Carac$dimension==0,]
  H1 <- Carac[Carac$dimension==1,]
  radios <- sort(unique(c(Carac$Birth,Carac$Death)))
  b0 <- c()
  b1 <- c()
   for (k in 1:length(radios)) {
     b0[k] <- nrow(H0[H0$Death > radios[k],])
     b1[k] <- nrow(H1[H1$Birth<=radios[k] & H1$Death>radios[k],])

   }
  Xcar <- b0-b1
  data.frame(radios,Xcar)
}
#X_carac(Y,X[,1],0.2)
#barcode_plotter(Y,X,0.2)
#(El barcode se hace con rips lo hace con Rips)



# Función de la caracteristica de Euler(todos los datos)-----

# Xdat<-Matriz de todos los datos X

Xcars <- function(Ydat,Xdat,maxscale=0.2,steps=pi/6) {
  pru <- list()
  lg <- list()


  pru <- lapply(
    X = seq_len(ncol(Xdat)),
    FUN = function(i) {
      Seccion_angulos(Y, X[, i], steps)
    }
  )#Data por angulos


  lg <- lapply(
    X = seq_len(ncol(Xdat)),
    FUN = function(j) {
      lapply(
        X = seq_len(length(pru[[1]])),
        FUN = function(i) {
          df <-
            X_carac(pru[[j]][[i]][, "y"], pru[[j]][[i]][, "x"], maxscale)
          cbind(df, angle = pru[[j]][[i]][1, "angle"])
            }
      )

    }
  )
  lg#Lista grande

}

E_car <- Xcars(Ydat,Xdat, steps = pi/4)#Cada numero de sublista representa un angulo y cada lista X_i
# Gráficos----

df1 <-
  bind_cols(bind_rows(E_car[[3]], .id = "angulo"), side = "upper")

df2 <-
  bind_cols(bind_rows(rev(E_car[[3]]), .id = "angulo"), side = "lower")

df <- bind_rows(df1, df2)

df <-  df %>%
  mutate(angulo = fct_reorder(angulo, as.numeric(angulo)))


ggplot(df,
       aes(radios, Xcar, colour = side)) +
  geom_line(size = 2) +
  # scale_color_viridis_d() +
  facet_wrap(. ~ angulo) +
  xlim(c(0, 0.003)) +
  theme_minimal()

df <- dplyr::bind_rows(E_car[[1]], .id = "angulo")

df <-  df %>%
  mutate(angulo = fct_reorder(angulo, as.numeric(angulo)))

ggplot(df,
       aes(radios, Xcar, colour = angulo)) +
  geom_line() +
  scale_color_viridis_d()+
  xlim(c(0,0.01))

df <- dplyr::bind_rows(E_car[[3]], .id = "angulo")

df <-  df %>%
  mutate(angulo = fct_reorder(angulo, as.numeric(angulo)))

ggplot(df,
       aes(radios, Xcar, colour = angle )) +
  geom_point() +
  scale_color_viridis_c()+
  xlim(c(0,0.005))





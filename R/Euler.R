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
    #names(angulos_data) <- cita
    angulos_data

  }

# pru <- Seccion_angulos(Y, X[, 1], pi / 6)
# ggplot() + geom_point(aes(X[, 1] - mean(X[, 1]), Y - mean(Y)), shape = 2) + geom_point(aes(pru[[2]][, 1], pru[[2]][, 2]))
# length(pru)

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
pru[[1]][[1]][,"y"]
X_carac(pru[[1]][[8]][,"y"],pru[[1]][[8]][,"x"],0.2)
#barcode_plotter(Y,X,0.2)
#(El barcode se hace con rips lo hace con Rips)



# Función de la caracteristica de Euler(todos los datos)-----

# Xdat<-Matriz de todos los datos X

Xcars <- function(Ydat,Xdat,maxscale=0.2,steps=pi/6) {
  #Etiquetas
  vect_angulos<- seq(0,2*pi,pi/24)#Se hizo con 24
  vect_angulos_caracter<- c('0','pi/24','pi/12','pi/8','pi/6','5*pi/24','pi/4','7*pi/24','pi/3','3*pi/8','5*pi/12','11*pi/24','pi/2','13*pi/24','7*pi/12','5*pi/8','2*pi/3','17*pi/24','3*pi/4','19*pi/24','5*pi/6','7*pi/8','11*pi/12','23*pi/24','pi','25*pi/24','13*pi/12','9*pi/8','7*pi/6','29*pi/24','5*pi/4','31*pi/24','4*pi/3','11*pi/8','17*pi/12','35*pi/24','3*pi/2','37*pi/24','19*pi/12','13*pi/8','5*pi/3','41*pi/24','7*pi/4','43*pi/24','11*pi/6','15*pi/8','23*pi/12','47*pi/24','2*pi')
  cita <- seq(0,2*pi,steps)
  indice <- c()
  for (d in 1:length(cita)) {
    indice[d] <- which(abs(vect_angulos-cita[d])<0.0001)
  }
  cita_caracter <- vect_angulos_caracter[indice]
  #library(rSymPy)
  # cita <- seq(0,2*pi,steps)
  # b <- paste("nsimplify(" ,cita,",[pi], tolerance=0.0000000001)")
  # vector <- lapply(b, sympy)
  # cita_caracter <-unlist(vector)
  # Verficar si lo anterior corre
  upper<- cita_caracter[which((cita-pi)<0)]
  lower <- cita_caracter[which(cita-pi>=0)]
  lower <- lower[which(lower!="2*pi")]
  Etiquetas <- rep(paste(upper,lower,sep = "--"),2)
  #Finaliza etiquetas

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
          cbind(df, angle = pru[[j]][[i]][1, "angle"]) %>%mutate(angle_numeric=cita[i])%>%  mutate(angle=cita_caracter[i]) %>% mutate(side=ifelse(angle%in%upper,"upper","lower")) %>% mutate(Etiqueta=Etiquetas[i])
        }
      )

    }
  )
   # for (k in 1:ncol(Xdat)) {
   #   names(lg[[k]]) <- cita_caracter
   # }
  lg#Lista grande

}





E_car <- Xcars(Ydat,Xdat, steps = pi/6)#Cada numero de sublista representa un angulo y cada lista X_i

# Gráficos----
df <- bind_rows(E_car)
#
# df1 <-
#   bind_cols(bind_rows(E_car[[3]], .id = "angulo"), side = "upper")
#
# df2 <-
#   bind_cols(bind_rows(rev(E_car[[3]]), .id = "angulo"), side = "lower")
#
# df <- bind_rows(df1, df2)
#
# df <-  df %>%
#   mutate(angulo = fct_reorder(angulo, as.numeric(angulo)))
# df$angulo
# unique(df$angle)
#
 ggplot(df,
        aes(radios, Xcar, colour = side)) +
   geom_line(size = 2) +
   #scale_color_viridis_d() +
   facet_wrap(. ~ Etiqueta) +
   xlim(c(0, 0.03)) +
  theme_minimal()


# df <- dplyr::bind_rows(E_car[[1]], .id = "angulo")
#
 df <-  df %>%
   mutate(angle = fct_reorder(angle, angle_numeric))

ggplot(df,
       aes(radios, Xcar, colour = angle)) +
  geom_line() +
  scale_color_viridis_d()+
  xlim(c(0,0.01))

# df <- dplyr::bind_rows(E_car[[3]], .id = "angulo")
#
# df <-  df %>%
#   mutate(angulo = fct_reorder(angulo, as.numeric(angulo)))

ggplot(df,
       aes(radios, Xcar, colour = angle)) +
  geom_point() +
  scale_color_viridis_c()+
  xlim(c(0,0.005))

############

# #Funcion para seccionar los datos por niveles ----
# #
# Seccion_niveles <- function(Yr,Xr,steps=10){
#   X_i <- Xr - mean(Xr)
#   Y_i <- Yr - mean(Yr)
#   Medida <- abs(max(Y_i) -min(Y_i))/steps
#   niveles <- seq(min(Y_i),max(Y_i),by=Medida)
#   niveles_data <- list(NULL)
#   niveles_data <- lapply(X=1:length(niveles),function(k){cbind(x=X_i[which(Y_i>=niveles[k])],y=Y_i[which(Y_i>=niveles[k])],nivel=niveles[k])})
#   niveles_data
# }
#
# #pru <- Seccion_niveles(Y,X[,1],10)
# # library(ggplot2)
# # ggplot()+geom_point(aes(X[,1] - mean(X[,1]),Y-mean(Y)),shape=2)+ geom_point(aes(pru[[2]][,1],pru[[2]][,2]))
#
#
#
# # Función de la característica de Euler con niveles----
#
# Xcars_niveles <- function(Ydat,Xdat,maxscale=0.2,steps=10) {
#   pru <- lapply(
#     X = seq_len(ncol(Xdat)),
#     FUN = function(i) {
#       Seccion_niveles(Y, X[,i], steps)
#     }
#   )#Data por angulos
#
#   lg <- lapply(
#     X = seq_len(ncol(Xdat)),
#     FUN = function(j) {
#       lapply(
#         X = seq_len(length(pru[[1]])),
#         FUN = function(i) {
#             X_carac(pru[[j]][[i]][, "y"], pru[[j]][[i]][, "x"], maxscale)
#
#         }
#       )
#
#     }
#   )
#   # for (k in 1:ncol(Xdat)) {
#   #   names(lg[[k]]) <- cita_caracter
#   # }
#   lg#Lista grande
#
# }
#
# E_carniveles <- Xcars_niveles(Ydat,Xdat, steps =10)

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
Seccion_angulos <- function(Yr,Xr,steps=pi/6){ # Reciben vector Y, X_i y el paso
      X_i <- Xr - mean(Xr)
      Y_i <- Yr - mean(Yr)
      angulos <- atan2(Y_i,X_i)
      angulos <- ifelse(angulos>=0,angulos,angulos + 2*pi)
      cita <- seq(0, 2*pi - steps, by=steps)
      angulos_data <- list(NULL)
      for (l in 1:length(cita)) {
        if(cita[l]<=pi){
          angulos_data[[l]] <- cbind(X_i[which(angulos>cita[l] & angulos<cita[l]+pi)],Y_i[which(angulos>cita[l] & angulos<cita[l]+pi)])
        }else{
          angulos_data[[l]] <- cbind(X_i[which(angulos>=cita[l] | angulos<=-pi+cita[l])],Y_i[which(angulos>=cita[l] | angulos<=-pi+cita[l])])
        }


      }
      names(angulos_data) <- paste("cita",0:(length(cita)-1),sep = "")
      angulos_data
}

pru <- Seccion_angulos(Y,X[,1],pi/2)
ggplot()+geom_point(aes(X[,1] - mean(X[,1]),Y-mean(Y)),shape=2)+ geom_point(aes(pru[[1]][,1],pru[[1]][,2]))


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

    ord <-
      order(d$diagram[, "Death"] - d$diagram[, "Birth"], decreasing = TRUE)

    dim1 <- d$diagram[ord, "dimension"] == 1
    d$diagram[head(ord[dim1]), "dimension"] <- -1
d$diagram
}
#caracteristica(Y,X[,1],0.2)


#Esta función recibe la columna  Y, X_i
X_carac <- function(Y,X,maxscale) {
  Carac <- as.data.frame(caracteristica(Ydat = Y, Xdat = X, maxscale))
  Carac <- Carac[Carac$dimension!=-1,]
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


  pru <- lapply(X = seq_len(ncol(Xdat)),
        FUN = function(i) {
         Seccion_angulos(Y,X[,i],steps)})#Data por angulos


  lg <- lapply(X = seq_len(ncol(Xdat)),
               FUN = function(j){
    lapply(X = seq_len(length(pru[[1]])),
              FUN = function(i) {
                X_carac(pru[[j]][[i]][,2],pru[[j]][[i]][,1],maxscale)})
               })
  lg#Lista grande

}

E_car <- Xcars(Ydat,Xdat)#Cada numero de sublista representa un angulo y cada lista X_i
# Gráficos----
ggplot(dplyr::bind_rows(E_car[[1]], .id="angulo"), aes(radios, Xcar, colour=angulo)) +
  geom_line()
ggplot(dplyr::bind_rows(E_car[[2]], .id="angulo"), aes(radios, Xcar, colour=angulo)) +
  geom_line()
ggplot(dplyr::bind_rows(E_car[[3]], .id="angulo"), aes(radios, Xcar, colour=angulo)) +
  geom_line()




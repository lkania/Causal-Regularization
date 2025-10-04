###########################################################################
# Compute moments
###########################################################################

moments <- function(data) {

  ne <- dim(data$Xe)[1]
  XYe <- t(data$Xe) %*% data$ye / ne
  XXe <- t(data$Xe) %*% data$Xe / ne

  no <- dim(data$Xo)[1]
  XYo <- t(data$Xo) %*% data$yo / no
  XXo <- t(data$Xo) %*% data$Xo / no

  Zplus <- XYe + XYo
  Gplus <- XXe + XXo

  Z <- XYe - XYo
  G <- XXe - XXo

  GZ <- t(G) %*% Z
  GG <- t(G) %*% G

  return(list(Xe = data$Xe,
              ye = data$ye,
              Xo = data$Xo,
              yo = data$yo,
              XYe = XYe,
              XXe = XXe,
              ne = ne,
              XYo = XYo,
              XXo = XXo,
              no = no,
              Z = Z,
              GZ = GZ,
              G = G,
              GG = GG,
              Gplus = Gplus,
              Zplus = Zplus,
              p = dim(data$Xe)[2]))
}
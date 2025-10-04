####################################################################
# Your R session's cwd should be the root of the project
####################################################################
####################################################################
# Causal Regularization with absolute l2-norm regularizer
####################################################################

source("./src/datautil.R") # load data utilities


compute_scr <- function(m, l) {

  # if (l == 1) {
  #   Gl <- m$G
  #   Zl <- m$Z
  # } else {
  Gl <- (1 - l) * m$Gplus + l * m$GG
  Zl <- (1 - l) * m$Zplus + l * m$GZ
  # }

  # Solve uses the Moore-Penrose generalized inverse of matrix A
  # (function ginv from package MASS).
  # See https://cran.r-project.org/web/packages/limSolve/limSolve.pdf
  return(limSolve::Solve(A = Gl, B = Zl, tol = 1e-10)) #1e-15

}

scr <- function(data, length, gamma = NULL) {

  m <- moments(data)

  if (is.null(gamma)) {
    gamma <- seq(0, 1, length = length)
  }

  beta <- NULL

  for (l in gamma) {

    beta <- rbind(beta, t(compute_scr(m = m, l = l)))

  }

  return(list(gamma = gamma, beta = beta, name = "CR_L2"))

}
###########################################################################
# Data generation
###########################################################################

data_gen <- function(target, B, n,
                     obs_mean = 0,
                     intervention_mean = 1,
                     shift_strength = 0,
                     shift_target = 0,
                     confounding = TRUE,
                     add_intercept = TRUE) {

  do <- gen(target = target,
            B = B,
            n = n,
            mean = obs_mean,
            shift_strength = 0,
            confounding = confounding,
            shift_target = 0,
            add_intercept = add_intercept)

  de <- gen(target = target,
            B = B,
            n = n,
            mean = intervention_mean,
            shift_strength = shift_strength,
            confounding = confounding,
            shift_target = shift_target,
            add_intercept = add_intercept)

  return(list(Xe = de$X, ye = de$y, Xo = do$X, yo = do$y))

}

gen <- function(target,
                B,
                n,
                mean = 0,
                shift_strength = 0,
                shift_target = 0,
                confounding = TRUE,
                tol = 1e-15,
                add_intercept = TRUE) {

  p <- dim(B)[1]

  distribution_shift <- 0
  if (shift_strength > 0) {

    distribution_shift <- sqrt(shift_strength) * matrix(
      rnorm(p * n, mean = mean, sd = 1),
      nrow = p)

    if (shift_target > 0) {
      distribution_shift[target,] <- sqrt(shift_target) * distribution_shift[target,]
    }else {
      distribution_shift[target,] <- 0
    }

  }

  confounding_ <- 0
  if (confounding) {
    confounding_ <- t(matrix(1, n, p) * rnorm(n, mean = 0, sd = 1))
  }

  noise <- matrix(rnorm(p * n, mean = 0, sd = 1), nrow = p)

  all <- t(limSolve::Solve(A = diag(p) - B,
                           B = distribution_shift + noise + confounding_,
                           tol = tol))

  X <- all[, -target]
  y <- as.matrix(all[, target])

  if (add_intercept) {
    X <- cbind(1, X)
  }


  return(list(X = X, y = y))
}
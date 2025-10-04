###########################################################################
# Measurements
###########################################################################

risk_vector <- function(X, y, beta) {
  return((y - X %*% beta)^2)
}

risk <- function(X, y, beta) {
  n <- dim(X)[1]
  pred <- as.vector(X %*% beta)
  y <- as.vector(y)
  return(sum((y - pred)^2) / n)
}

risk_on_data <- function(data, beta) {
  return(risk(X = data$X, y = data$y, beta = beta))
}

risk_on_intervention_data <- function(data, beta) {
  return(risk(X = data$Xe, y = data$ye, beta = beta))
}

risk_on_obs_data <- function(data, beta) {
  return(risk(X = data$Xo, y = data$yo, beta = beta))
}

in_sample_risk <- function(data, beta) {
  return(risk_on_intervention_data(data = data, beta = beta) + risk_on_obs_data(data = data, beta = beta))
}

risk_difference <- function(data, beta) {
  return(risk_on_intervention_data(data = data, beta = beta) - risk_on_obs_data(data = data, beta = beta))
}

absolute_risk_difference <- function(data, beta) {
  return(abs(risk_difference(data, beta)))
}

squared_risk_difference <- function(data, beta) {
  return((risk_difference(data, beta))^2)
}

norm2 <- function(beta) {
  return(sqrt(sum(beta * beta)))
}

measure <- function(data, beta, metric) {

  n <- dim(beta)[1]
  measurements <- numeric(n)

  for (i in seq(1, n, 1)) {

    b <- t(t(beta[i,]))

    measurements[[i]] <- metric(data = data, beta = b)
  }

  return(measurements)

}
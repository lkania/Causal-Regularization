source('./src/measurements.R')

###########################################################################
# Fit estimators using all data / sample splitting / cross-validation
###########################################################################

# Note: The following function assumes that the range of gamma is fixed
fit_cross_validation <- function(folds,
                                 f,
                                 data_train,
                                 data_val,
                                 metric,
                                 plot = FALSE) {

  risk_per_fold <- vector("list", folds)

  for (i in 1:folds) {

    fsol <- f(data = data_train[[i]])
    gammas <- fsol$gamma

    risk_per_fold[[i]] <- measure(
      beta = fsol$beta,
      data = data_val[[i]],
      metric = metric)

  }

  risks <- do.call(rbind, risk_per_fold)

  risk_mean <- apply(risks, 2, mean)
  risk_sd <- apply(risks, 2, sd)

  stopifnot(length(risk_mean) == length(gammas))
  stopifnot(length(risk_sd) == length(gammas))

  idx_selected <- which.min(risk_mean)
  gamma_selected <- gammas[idx_selected]

  ########################################
  # plot cross-validation selection
  ########################################
  p <- NULL
  if (plot) {
    df <- data.frame(x = gammas, y = risk_mean,
                     fold = 0, type = 'Average')
    for (i in 1:folds) {
      df <- rbind(df, data.frame(
        x = gammas,
        y = risk_per_fold[[i]],
        fold = i, type = 'Fold'))
    }
    df <- rbind(df, data.frame(
      x = c(gamma_selected, gamma_selected),
      y = c(min(df$y), max(df$y)),
      fold = -1, type = 'Selected'))

    p <- ggplot(df) +
      geom_line(mapping =
                  aes(x = x, y = y, group = fold, color = type, linetype = type)) +
      scale_color_manual(values = c('black', 'black', 'red')) +
      labs(color = "", linetype = "") +
      theme(legend.title = element_blank())

    ticks <- c(0, 0.25, 0.5, 0.75, 1, gamma_selected)
    label_ticks <- expression(
      OLS, '0.25', '0.5', '0.75', CD, CR)
    idx <- order(ticks)

    p <- p +
      scale_x_continuous(breaks = ticks[idx], labels = label_ticks[idx]) +
      ylab('Absolute risk difference') +
      xlab(expression('Regularization path' ~ tilde(beta)[lambda]))
  }


  return(list(opt = gamma_selected,
              risk_mean = risk_mean,
              risk_sd = risk_sd,
              plot = p))
}

subset_data <- function(data, sele = NULL, selo = NULL) {

  Xo <- data$Xo
  yo <- data$yo

  Xe <- data$Xe
  ye <- data$ye

  if (!is.null(selo)) {
    Xo <- as.matrix(Xo[selo,])
    yo <- as.vector(yo[selo])
  }

  if (!is.null(sele)) {
    Xe <- as.matrix(Xe[sele,])
    ye <- as.vector(ye[sele])
  }

  return(list(Xe = Xe, ye = ye, Xo = Xo, yo = yo))

}

split_idx <- function(folds, n) {

  #equal sized splits
  g <- sample(rep(1:folds, floor(n / folds)))

  res <- n %% folds
  if (res > 0) {
    g <- c(g, sample(1:folds, res))
  }

  return(g)

}

cross_validation <- function(folds, data, estimators, metric, plot = FALSE) {
  # The functions assumes that all provided estimators have fixed range

  go <- split_idx(folds = folds, n = dim(data$Xo)[1])
  ge <- split_idx(folds = folds, n = dim(data$Xe)[1])

  data_train <- vector("list", folds)
  data_val <- vector("list", folds)

  for (i in 1:folds) {

    data_train[[i]] <- subset_data(data,
                                   sele = (ge != i),
                                   selo = (go != i))

    data_val[[i]] <- subset_data(data,
                                 sele = (ge == i),
                                 selo = (go == i))

  }

  estimands <- vector("list", length(estimators))

  for (k in 1:length(estimators)) {

    estimands[[k]] <- fit_cross_validation(folds = folds,
                                           metric = metric,
                                           f = estimators[[k]],
                                           data_train = data_train,
                                           data_val = data_val,
                                           plot = plot)

  }

  return(estimands)
}
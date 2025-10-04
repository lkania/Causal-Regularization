# assuming it's run from the experiment folder
source("./src/selection.R")
source("./src/measurements.R") # load metrics
source("./src/plotutils.R")

library(dplyr)
library(tidyr)

select <- function(data, idx, names) {
  data_ <- data[idx,]

  X <- as.matrix(data_[, (colnames(data_) %in% names)])
  y <- as.vector(data_$y)

  stopifnot(length(y) >= 1)

  return(list(X = X, y = y))
}

center_train_data <- function(data) {

  data_ <- data
  meanY <- mean(data_$yo)
  meanX <- mean(data_$Xo)

  data_$ye <- data_$ye - meanY
  data_$yo <- data_$yo - meanY
  data_$Xe <- data_$Xe - meanX
  data_$Xo <- data_$Xo - meanX

  return(list(data = data_, meanY = meanY, meanX = meanX))

}

resample <- function(data_train,
                     data_test,
                     compute_estimator_for_all_gamma,
                     folds_for_cv,
                     metric_for_cv,
                     path_to_save_results = NULL,
                     n_resamples = NULL,
                     prefix = NULL) {

  if (is.null(n_resamples)) {
    n_resamples <- 1
  }

  if (is.null(prefix)) {
    prefix <- ""
  }else {
    prefix <- paste0(prefix, "_")
  }

  centering_data_train <- center_train_data(data_train)
  centred_data_train <- centering_data_train$data
  centred_data_test_based_on_data_train <- list(
    X = data_test$X - centering_data_train$meanX,
    y = data_test$y - centering_data_train$meanY)

  beta_path_on_data_train <- compute_estimator_for_all_gamma(
    data = centred_data_train)

  method_name <- beta_path_on_data_train$name
  suffix <- method_name

  out_of_sample_risk_on_data_test <- measure(
    data = centred_data_test_based_on_data_train,
    beta = beta_path_on_data_train$beta,
    metric = risk_on_data)

  gamma_opt_on_data_test <- beta_path_on_data_train$gamma[which.min(out_of_sample_risk_on_data_test)]

  stopifnot(length(out_of_sample_risk_on_data_test) == length(beta_path_on_data_train$gamma))

  run_resample <- function(resample_idx) {


    centred_resampled_data_train <- centred_data_train
    centering_resampled_data_train <- centering_data_train

    if (n_resamples > 1) {

      data_intervention_train_idx <- sample(1:length(data_train$ye), replace = FALSE)
      data_obs_train_idx <- sample(1:length(data_train$yo), replace = FALSE)

      resampled_data_train <- subset_data(data = data_train,
                                          sele = data_intervention_train_idx,
                                          selo = data_obs_train_idx)
      centering_resampled_data_train <- center_train_data(resampled_data_train)
      centred_resampled_data_train <- centering_resampled_data_train$data
    }

    estimands <- cross_validation(
      folds = folds_for_cv,
      data = centred_resampled_data_train,
      estimators = c(compute_estimator_for_all_gamma),
      metric = metric_for_cv,
      plot = FALSE)

    in_sample_risk_mean <- estimands[[1]]$risk_mean
    in_sample_risk_sd <- estimands[[1]]$risk_sd
    gamma_star <- estimands[[1]]$opt

    stopifnot(length(in_sample_risk_sd) == length(beta_path_on_data_train$gamma))
    stopifnot(length(in_sample_risk_mean) == length(beta_path_on_data_train$gamma))

    if (n_resamples > 1) {
      beta_path_on_resampled_data_train <- compute_estimator_for_all_gamma(
        data = centred_resampled_data_train)

      centred_data_test_based_on_resampled_data_train <- list(
        X = data_test$X - centering_resampled_data_train$meanX,
        y = data_test$y - centering_resampled_data_train$meanY)

      out_of_sample_risk <- measure(
        data = centred_data_test_based_on_resampled_data_train,
        beta = beta_path_on_resampled_data_train$beta,
        metric = risk_on_data)

    }else {

      beta_path_on_resampled_data_train <- beta_path_on_data_train
      out_of_sample_risk <- out_of_sample_risk_on_data_test

    }

    df <- data.frame(gamma_star = gamma_star,
                     in_sample_risk_mean = in_sample_risk_mean,
                     in_sample_risk_sd = in_sample_risk_sd,
                     out_of_sample_risk = out_of_sample_risk,
                     gammas = beta_path_on_resampled_data_train$gamma,
                     idx = resample_idx)
    df$idx <- as.factor(df$idx)

    return(df)

  }

  if (n_resamples > 1) {

    # The following line is a parallel version of
    # results <- lapply(X = 1:n_resamples, FUN = run)
    # results <- future_lapply(X = 1:n_resamples, FUN = run, future.seed = TRUE)
    # df <- do.call(rbind, results)
    df <- future_lapply(1:n_resamples, run_resample, future.seed = TRUE) %>% bind_rows()

    # Furthermore, note that future_lapply protects against recurvise parallelism
    # See: "Built-in protection against recursive parallelism" in
    # https://cran.r-project.org/web/packages/future/vignettes/future-3-topologies.html

    # summary of out-of-sample risk
    means <- aggregate(out_of_sample_risk ~ gammas, df, mean)
    colnames(means)[colnames(means) == "out_of_sample_risk"] <- "out_of_sample_risk_mean"
    sds <- aggregate(out_of_sample_risk ~ gammas, df, sd)
    colnames(sds)[colnames(sds) == "out_of_sample_risk"] <- "out_of_sample_risk_sd"
    out_of_sample_risk_summary <- merge(means, sds, by = "gammas")
    df <- merge(df, out_of_sample_risk_summary, by = "gammas")

    # summary of model selection risk
    means <- aggregate(in_sample_risk_mean ~ gammas, df, mean)
    colnames(means)[colnames(means) == "in_sample_risk_mean"] <- "in_sample_risk_mean"
    sds <- aggregate(in_sample_risk_mean ~ gammas, df, sd)
    colnames(sds)[colnames(sds) == "in_sample_risk_mean"] <- "in_sample_risk_sd"
    in_sample_risk_summary <- merge(means, sds, by = "gammas")
    df <- df[, !(colnames(df) %in% c(
      "in_sample_risk_mean",
      "in_sample_risk_sd"))]
    df <- merge(df, in_sample_risk_summary, by = "gammas")


  }else {

    # No resample is done
    df <- run_resample(resample_idx = 1)
    df$out_of_sample_risk_mean <- df$out_of_sample_risk
    df$out_of_sample_risk_sd <- 0

  }

  df <- df %>% rename(gamma_selected = gamma_star)
  df$gamma_opt <- gamma_opt_on_data_test


  if (!is.null(path_to_save_results)) {

    writeLines(paste0("Train observation dataset n=", length(data_train$yo)))
    writeLines(paste0("Train intervention dataset n=", length(data_train$ye)))
    writeLines(paste0("Test dataset n=", length(data_test$y)))

    gamma_stats <- df %>%
      # there are many rows per idx, all of them have the same
      # gamma_selected and gamma_opt
      group_by(idx) %>%
      summarise(
        gamma_cr = mean(gamma_selected),
        gamma_opt = mean(gamma_opt),
        .groups = "drop"
      ) %>%
      pivot_longer(cols = c(gamma_cr, gamma_opt),
                   names_to = "label",
                   values_to = "gamma") %>%
      mutate(label = recode(label,
                            gamma_cr = "CR",
                            gamma_opt = "OPT")) %>%
      group_by(label) %>%
      summarise(
        gamma_mean = mean(gamma),
        gamma_sd = sd(gamma),
        gamma_lb = pmax(gamma_mean - gamma_sd, 0),
        gamma_ub = pmin(gamma_mean + gamma_sd, 1),
        .groups = "drop"
      )

    #
    # Plot the model selection loss and the out-of-sample risk together
    #

    df_long <- df %>%
      dplyr::select(gammas,
                    in_sample_risk_mean, in_sample_risk_sd,
                    out_of_sample_risk_mean, out_of_sample_risk_sd) %>%
      pivot_longer(
        cols = c(in_sample_risk_mean, out_of_sample_risk_mean,
                 in_sample_risk_sd, out_of_sample_risk_sd),
        names_to = c("type", ".value"),
        names_pattern = "(in_sample|out_of_sample)_risk_(.*)"
      ) %>%
      mutate(
        type = recode(type,
                      in_sample = "Cross-validated in-sample risk",
                      out_of_sample = "Out-of-sample risk")
      )

    p <- ggplot(df_long, aes(x = gammas, y = mean)) +
      geom_line(color = 'black') +
      geom_ribbon(mapping = aes(ymin = mean - sd, ymax = mean + sd),
                  alpha = 0.1,
                  linewidth = 0,
                  fill = 'black') +
      geom_vline(data = gamma_stats,
                 mapping = aes(xintercept = gamma_mean, color = label),
                 linetype = "dashed",
                 size = 0.7) +
      geom_rect(data = gamma_stats,
                mapping = aes(xmin = gamma_lb,
                              xmax = gamma_ub,
                              ymin = -Inf, ymax = Inf,
                              fill = label),
                alpha = 0.1,
                inherit.aes = FALSE) +
      facet_wrap(~type, scales = "free_y") +
      scale_x_continuous(
        breaks = seq(0, 1, 0.25),
        labels = expression(OLS, '0.25', '0.5', '0.75', CD)
      ) +
      scale_fill_discrete(name = "Method") +
      scale_color_discrete(name = "Method") +
      ylab(NULL) +
      xlab(expression('Causal regularization path' ~ tilde(beta)[lambda])) +
      theme(legend.position = "top",
            strip.text.x = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12)) +
      guides(color = guide_legend(nrow = 1))

    save_plot(p,
              paste0(path_to_save_results, prefix, "result_", suffix, ".pdf"),
              dims = list(width = 10, height = 4))

  }

  return(df)

}
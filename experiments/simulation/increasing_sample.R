####################################################################
# Your R session's working directory should be the root of the project
####################################################################
start_time <- Sys.time()

# Read arguments
source("./src/parse.R")
batch_start <- as.integer(get_arg_value("batchStart", NA))
batch_end <- as.integer(get_arg_value("batchEnd", NA))
shift_target <- as.double(get_arg_value("shiftTarget", 0))
in_sample_shift <- as.double(get_arg_value("inSampleShift", 5))
intervention_mean <- as.double(get_arg_value("interventionMean", 0))

source("./src/init_parallel.R")
source("./src/data.R")
source("./src/utils.R")
source("./src/measurements.R")
source("./src/selection.R")
source("./src/real.R")
source("./experiments/simulation/simulation.R")
source("./src/cr_l2.R")
path <- './experiments/simulation/'

library(digest)
library(purrr)

# Cache parameters
batch_size <- 100

# Parameters used by cross-validation
metric_for_cv <- absolute_risk_difference
compute_estimator_for_all_gamma <- (function(data) scr(data = data, length = 50))
folds_for_cv <- 10

# Simulation parameters
out_of_sample_size <- 1e6
in_sample_sizes <- c(1e2, 1e3, 1e4)
in_sample_shifts <- c(in_sample_shift)
out_of_sample_shifts <- c(10, 50, 100)
in_sample_obs_mean <- 0
confounding <- TRUE
target <- 4
B <- matrix(c(0, 0, 0, 0, 0, 0, 0,
              1, 0, 0, 0, 0, 0, 0,
              1, 1, 0, 0, 0, 0, 0,
              0, 1, 1, 0, 0, 0, 0,
              0, 0, 0, 1, 0, 0, 0,
              0, 0, 0, 1, 0, 0, 0,
              0, 0, 0, 0, 0, 1, 0),
            nrow = 7, byrow = TRUE)
n_trials <- 1000

# Directories
path <- paste0(path,
               "shift_target=",
               shift_target,
               "_intervention_mean=",
               intervention_mean,
               '_in_sample_shift=',
               in_sample_shifts,
               '_out_shift=',
               paste(out_of_sample_shifts, collapse = ""),
               "_")
checkpoint_dir <- paste0(path, "sample_checkpoints")
dir.create(checkpoint_dir, showWarnings = FALSE)


####################################################################
# Simulation
####################################################################

grid <- expand.grid(trial = 1:n_trials,
                    out_of_sample_shift = out_of_sample_shifts,
                    in_sample_shift = in_sample_shifts,
                    in_sample_size = in_sample_sizes)
stopifnot(!any(duplicated(grid)))

# Prepare batches
total_tasks <- nrow(grid)
batches <- split(grid, ceiling(seq_len(total_tasks) / batch_size))
total_batches <- length(batches)

# Find existing checkpoints
writeLines("Finding existing checkpoints.\n")
existing_files <- list.files(checkpoint_dir,
                             pattern = "^results_batch_\\d+\\.rds$")
done_batches <- as.integer(
  gsub("\\D",
       "",
       sub("results_batch_", "", existing_files)))

writeLines("Starting simulations.\n")
writeLines(paste0("Processing ", dim(grid)[1], " experiments.\n"))

overall_start_time <- Sys.time()

for (i in seq_along(batches)) {

  if ((!is.na(batch_start) && i < batch_start) || (!is.na(batch_end) && i > batch_end)) {
    next
  }

  checkpoint_file <- file.path(checkpoint_dir, paste0("results_batch_", i, ".rds"))

  if (file.exists(checkpoint_file)) {
    writeLines(sprintf(
      "Batch %d of %d already completed. Skipping.\n", i, total_batches))
    next
  }else {
    writeLines(sprintf("Processing Batch %d of %d\n", i, total_batches))
  }

  batch <- batches[[i]]
  batch_start_time <- Sys.time()

  # We set a unique seed per batch to avoid the order in which the batches
  # are processed affecting the final result of the simulation
  set.seed(i)

  # Use  future.seed = TRUE to produce sound random numbers
  # regardless of computing lapply sequentially or in parallel
  # See https://cran.r-project.org/web/packages/future.apply/future.apply.pdf
  result_batch <- future_mapply(
    FUN = function(...) {
      run_simulation(...) # must return a data.frame
    },
    trial = batch$trial,
    out_of_sample_shift = batch$out_of_sample_shift,
    in_sample_size = batch$in_sample_size,
    in_sample_shift = batch$in_sample_shift,
    MoreArgs = list(
      target = target,
      B = B,
      intervention_mean = intervention_mean,
      in_sample_obs_mean = in_sample_obs_mean,
      folds_for_cv = folds_for_cv,
      confounding = confounding,
      shift_target = shift_target
    ),
    SIMPLIFY = FALSE,
    future.seed = TRUE)

  # Merge all data frames from batch
  result_df <- do.call(rbind, result_batch)
  saveRDS(result_df, checkpoint_file)

  #Progress reporting
  batch_end_time <- Sys.time()
  duration <- as.numeric(difftime(batch_end_time, batch_start_time, units = "secs"))
  total_elapsed <- as.numeric(difftime(batch_end_time, overall_start_time, units = "secs"))
  avg_time <- total_elapsed / length(setdiff(seq_len(i), done_batches))
  batches_left <- total_batches - i
  eta_secs <- avg_time * batches_left
  eta_formatted <- format(.POSIXct(eta_secs, tz = "UTC"), "%H:%M:%S")

  writeLines(sprintf("Batch %d of %d completed (%.1f sec) | ETA: %s\n",
                     i, total_batches, duration, eta_formatted))
}

####################################################################
# Final merge of all batches
####################################################################
batch_files <- list.files(checkpoint_dir,
                          pattern = "^results_batch_\\d+\\.rds$",
                          full.names = TRUE)

if (length(batch_files) < total_batches) {

  writeLines(sprintf(
    "Only %d of %d batches completed. Exiting without merging.\n",
    length(batch_files), total_batches))

}else {

  writeLines("All batches completed. Proceeding to merge.\n")
  df <- do.call(rbind, lapply(batch_files, readRDS))
  df <- as.data.frame(df)

  df$out_of_sample_shift <- as.factor(df$out_of_sample_shift)
  df$in_sample_shift <- as.factor(df$in_sample_shift)
  df$in_sample_size <- as.factor(df$in_sample_size)
  df$out_of_sample_size <- as.factor(df$out_of_sample_size)


  ####################################################################
  # Pre-process data for plots
  ####################################################################

  gamma_stats <- df %>%
    group_by(in_sample_size, in_sample_shift) %>%
    summarise(
      gamma_mean = mean(gamma_selected),
      gamma_sd = sd(gamma_selected),
      gamma_lb = pmax(gamma_mean - gamma_sd, 0),
      gamma_ub = pmin(gamma_mean + gamma_sd, 1),
      .groups = "drop"
    )

  in_sample_risk_summary <- df %>%
    group_by(gammas, in_sample_size, in_sample_shift) %>%
    summarise(
      in_sample_risk_sd = sd(in_sample_risk_mean),
      in_sample_risk_mean = mean(in_sample_risk_mean),
      in_sample_risk_lb = pmax(in_sample_risk_mean - in_sample_risk_sd, 0),
      in_sample_risk_ub = in_sample_risk_mean + in_sample_risk_sd,
      .groups = 'drop'
    )

  out_of_sample_risk_summary <- df %>%
    group_by(gammas, in_sample_size, in_sample_shift, out_of_sample_shift) %>%
    summarise(
      out_of_sample_risk_mean = mean(out_of_sample_risk),
      out_of_sample_risk_sd = sd(out_of_sample_risk),
      out_of_sample_risk_lb = pmax(out_of_sample_risk_mean - out_of_sample_risk_sd, 0),
      out_of_sample_risk_ub = out_of_sample_risk_mean + out_of_sample_risk_sd,
      .groups = 'drop'
    )

  out_of_sample_risk_per_method <- function(gamma) {
    return(df %>%
             filter(gammas == gamma) %>%
             group_by(in_sample_size, in_sample_shift, out_of_sample_shift) %>%
             summarise(
               out_of_sample_risk_mean = mean(out_of_sample_risk),
               out_of_sample_risk_sd = sd(out_of_sample_risk),
               out_of_sample_risk_lb = pmax(out_of_sample_risk_mean - out_of_sample_risk_sd, 0),
               out_of_sample_risk_ub = out_of_sample_risk_mean + out_of_sample_risk_sd,
               .groups = 'drop'
             ))
  }


  method_summary <- imap_dfr(
    list(CR = df$gamma_selected, OLS = 0, CD = 1),
    function(gamma, method_name) {
      out_of_sample_risk_per_method(gamma) %>% mutate(method = method_name)
    })
  method_summary$method <- as.factor(method_summary$method)
  method_summary$in_sample_size <- as.numeric(as.character(method_summary$in_sample_size))
  method_summary$out_of_sample_shift <- factor(method_summary$out_of_sample_shift, levels = out_of_sample_shifts)

  ####################################################################
  # Plot performance of OLS, CR, and CD
  ####################################################################

  p1 <- ggplot(data = method_summary,
               mapping = aes(x = in_sample_size,
                             y = out_of_sample_risk_mean,
                             fill = method,
                             color = method)) +
    geom_line(alpha = 1) +
    geom_ribbon(mapping = aes(ymin = out_of_sample_risk_lb,
                              ymax = out_of_sample_risk_ub),
                alpha = 0.1,
                linewidth = 0) +
    facet_wrap(. ~ out_of_sample_shift,
               ncol = length(out_of_sample_shifts),
               scales = 'fixed',
               labeller = labeller(
                 out_of_sample_shift = function(labels) {
                   paste("Out-of-sample distribution shift =", labels)
                 }
               )) +
    scale_fill_discrete(name = "Method") +
    scale_color_discrete(name = "Method") +
    labs(
      x = 'In-sample size',
      y = 'Normalized out-of-sample risk'
    ) +
    theme(legend.position = "top",
          strip.text.x = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12)) +
    scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    coord_cartesian(ylim = c(0, 0.4))

  save_plot(p1,
            paste0(path, "increasing_sample_comparison.pdf"),
            dims = list(width = 10, height = 4))


  ####################################################################
  # Plot model selection loss
  ####################################################################

  p1 <- ggplot(data = in_sample_risk_summary,
               mapping = aes(x = gammas,
                             y = in_sample_risk_mean,
                             color = in_sample_shift,
                             fill = in_sample_shift)) +
    geom_line(alpha = 1) +
    geom_ribbon(mapping = aes(ymin = in_sample_risk_lb,
                              ymax = in_sample_risk_ub),
                alpha = 0.1,
                linewidth = 0) +
    facet_wrap(. ~ in_sample_size,
               ncol = length(in_sample_sizes),
               scales = 'fixed',
               labeller = labeller(
                 in_sample_size = function(labels) {
                   paste("In-sample size =", labels)
                 }
               )) +
    geom_vline(data = gamma_stats,
               aes(xintercept = gamma_mean),
               linetype = "dashed",
               size = 0.7) +
    geom_rect(data = gamma_stats,
              aes(xmin = gamma_lb,
                  xmax = gamma_ub,
                  ymin = -Inf, ymax = Inf),
              alpha = 0.1,
              inherit.aes = FALSE) +
    scale_x_continuous(
      breaks = seq(0, 1, 0.25),
      labels = expression(OLS, '0.25', '0.5', '0.75', CD)
    ) +
    ylab('Normalized in-sample absolute risk difference') +
    xlab(expression('Causal regularization path' ~ tilde(beta)[lambda])) +
    theme(legend.position = "none") +
    theme(strip.text.x = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)) +
    coord_cartesian(ylim = c(0, 1))

  save_plot(p1,
            paste0(path, "increasing_sample_in_sample_risk.pdf"),
            dims = list(width = 10, height = 4))


  ####################################################################
  # Plot out-of-sample risk
  ####################################################################

  p1 <- ggplot(data = out_of_sample_risk_summary,
               mapping = aes(x = gammas,
                             y = out_of_sample_risk_mean,
                             fill = out_of_sample_shift,
                             color = out_of_sample_shift)) +
    geom_line(alpha = 1) +
    geom_ribbon(mapping = aes(ymin = out_of_sample_risk_lb,
                              ymax = out_of_sample_risk_ub),
                alpha = 0.1,
                linewidth = 0) +
    facet_wrap(. ~ in_sample_size,
               ncol = length(in_sample_sizes),
               scales = 'fixed',
               labeller = labeller(
                 in_sample_size = function(labels) {
                   paste("In-sample size =", labels)
                 }
               )
    ) +
    geom_vline(data = gamma_stats,
               aes(xintercept = gamma_mean),
               linetype = "dashed",
               size = 0.7) +
    geom_rect(data = gamma_stats,
              aes(xmin = gamma_lb,
                  xmax = gamma_ub,
                  ymin = -Inf, ymax = Inf),
              alpha = 0.1,
              inherit.aes = FALSE) +
    scale_x_continuous(
      breaks = seq(0, 1, 0.25),
      labels = expression(OLS, '0.25', '0.5', '0.75', CD)
    ) +
    scale_fill_discrete(name = "Out-of-sample distribution shift") +
    scale_color_discrete(name = "Out-of-sample distribution shift") +
    labs(
      y = 'Normalized out-of-sample risk',
      x = expression('Causal regularization path' ~ tilde(beta)[lambda])) +
    theme(legend.position = "top") +
    theme(strip.text.x = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)) +
    coord_cartesian(ylim = c(0, 0.4))

  save_plot(p1,
            paste0(path, "increasing_sample_out_of_sample_risk.pdf"),
            dims = list(width = 10, height = 4))
}

####################################################################
# Print runtime
####################################################################
end_time <- Sys.time()
total_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
hours <- total_seconds %/% 3600
minutes <- (total_seconds %% 3600) %/% 60
seconds <- round(total_seconds %% 60)
writeLines(sprintf("Total runtime: %02d:%02d:%02d (hh:mm:ss).\n",
                   hours,
                   minutes,
                   seconds))
####################################################################
# Your R session's working directory should be the root of the project
####################################################################
start_time <- Sys.time()

# Read arguments
source("./src/parse.R")
batch_start <- as.integer(get_arg_value("batchStart", NA))
batch_end <- as.integer(get_arg_value("batchEnd", NA))
shift_target <- as.double(get_arg_value("shiftTarget", 0))
intervention_mean <- as.double(get_arg_value("interventionMean", 0))
out_of_sample_shift <- as.double(get_arg_value("outOfSampleShift", 100))

source("./src/init_parallel.R")
source("./src/data.R")
source("./src/utils.R")
source("./src/measurements.R")
source("./src/selection.R")
source("./src/real.R")
source("./experiments/simulation/simulation.R")
source("./src/cr_l2.R")
path <- './experiments/simulation/'

library(purrr)
library(dplyr)

# Cache parameters
batch_size <- 40

# Parameters used by cross-validation
metric_for_cv <- absolute_risk_difference
compute_estimator_for_all_gamma <- (function(data) scr(data = data, length = 50))
folds_for_cv <- 10

# Simulation parameters
out_of_sample_size <- 1e6
in_sample_shifts <- c(1)
out_of_sample_shifts <- c(out_of_sample_shift)
in_sample_obs_mean <- 0
confounding <- TRUE
in_sample_sizes <- c(1e2, 1e3, 1e4)
n_trials <- 30
dims <- c(10, 20, 50)
n_graphs <- 30

# Directories
path <- paste0(path,
               "shift_target=",
               shift_target,
               "_intervention_mean=",
               intervention_mean,
               '_out_shift=',
               out_of_sample_shifts,
               '_in_sample_sizes=',
               paste(in_sample_sizes, collapse = ""),
               '_dims=',
               paste(dims, collapse = ""),
               "_")
checkpoint_dir <- paste0(path, "dim_checkpoints")
dir.create(checkpoint_dir, showWarnings = FALSE)

# Every directed edge is sample with probability 1/2
uniform_random_graph <- function(dim, self_loops = FALSE) {
  B <- matrix(sample(0:1, dim * dim, replace = TRUE), nrow = dim, ncol = dim)
  if (!self_loops) {
    diag(B) <- 0
  }
  return(B)
}

####################################################################
# Simulation
####################################################################

grid <- expand.grid(
  trial = 1:n_trials,
  graph = 1:n_graphs,
  dim = dims,
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
existing_files <- list.files(checkpoint_dir, pattern = "^results_batch_\\d+\\.rds$")
done_batches <- as.integer(gsub("\\D", "", sub("results_batch_", "", existing_files)))


run_simulation_ <- function(trial,
                            graph,
                            dim,
                            out_of_sample_shift,
                            in_sample_size,
                            in_sample_shift
) {

  B <- uniform_random_graph(dim = dim)

  # choose node with highest in-degree as target
  target <- order(rowSums(B), decreasing = TRUE)[1]

  df <- run_simulation(trial = trial,
                       target = target,
                       B = B,
                       out_of_sample_shift = out_of_sample_shift,
                       in_sample_size = in_sample_size,
                       in_sample_shift = in_sample_shift,
                       intervention_mean = intervention_mean,
                       in_sample_obs_mean = in_sample_obs_mean,
                       folds_for_cv = folds_for_cv,
                       confounding = confounding,
                       shift_target = shift_target)

  df$dim <- dim
  df$graph <- graph
  return(df)

}

writeLines("Starting simulations.\n")
writeLines(paste0("Running ", dim(grid)[1], " experiments.\n"))

overall_start_time <- Sys.time()

for (i in seq_along(batches)) {

  if ((!is.na(batch_start) && i < batch_start) || (!is.na(batch_end) && i > batch_end)) {
    next
  }

  checkpoint_file <- file.path(checkpoint_dir, paste0("results_batch_", i, ".rds"))

  if (file.exists(checkpoint_file)) {
    writeLines(sprintf("Batch %d of %d already completed. Skipping.\n", i, total_batches))
    next
  }else {
    writeLines(sprintf("Processing Batch %d of %d\n", i, total_batches))
  }

  batch <- batches[[i]]
  batch_start_time <- Sys.time()

  # We set a unique seed per batch to avoid the order in which the batches
  # are processed affecting the final result of the simulation
  set.seed(i)

  result_batch <- future_mapply(
    FUN = function(...) {
      run_simulation_(...)
    },
    dim = batch$dim,
    graph = batch$graph,
    trial = batch$trial,
    out_of_sample_shift = batch$out_of_sample_shift,
    in_sample_size = batch$in_sample_size,
    in_sample_shift = batch$in_sample_shift,
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


# Final merge of all batches
batch_files <- list.files(checkpoint_dir,
                          pattern = "^results_batch_\\d+\\.rds$",
                          full.names = TRUE)

if (length(batch_files) < total_batches) {
  writeLines(sprintf(
    "Only %d of %d batches completed. Exiting without merging.\n",
    length(batch_files), total_batches))
} else {

  writeLines("All batches completed. Proceeding to merge.\n")

  df <- do.call(rbind, lapply(batch_files, readRDS))
  df <- as.data.frame(df)

  df$out_of_sample_shift <- as.factor(df$out_of_sample_shift)
  df$in_sample_shift <- as.factor(df$in_sample_shift)
  df$in_sample_size <- as.factor(df$in_sample_size)
  df$out_of_sample_size <- as.factor(df$out_of_sample_size)
  df$dim <- as.factor(df$dim)

  ####################################################################
  # Pre-process data for plots
  ####################################################################

  gamma_stats <- df %>%
    group_by(in_sample_size, in_sample_shift, dim) %>%
    summarise(
      gamma_mean = mean(gamma_selected),
      gamma_sd = sd(gamma_selected),
      .groups = "drop"
    )

  in_sample_risk_summary <- df %>%
    group_by(gammas, in_sample_size, in_sample_shift, dim) %>%
    summarise(
      in_sample_risk_sd = sd(in_sample_risk_mean),
      in_sample_risk_mean = mean(in_sample_risk_mean),
      .groups = 'drop'
    )

  out_of_sample_risk_summary <- df %>%
    group_by(gammas, in_sample_size, in_sample_shift, out_of_sample_shift, dim) %>%
    summarise(
      out_of_sample_risk_mean = mean(out_of_sample_risk),
      out_of_sample_risk_sd = sd(out_of_sample_risk),
      .groups = 'drop'
    )

  out_of_sample_risk_per_method <- function(gamma) {
    return(df %>%
             filter(gammas == gamma) %>%
             group_by(in_sample_size, in_sample_shift, out_of_sample_shift, dim) %>%
             summarise(
               out_of_sample_risk_mean = mean(out_of_sample_risk),
               out_of_sample_risk_sd = sd(out_of_sample_risk),
               .groups = "drop"
             ))
  }


  method_summary <- imap_dfr(
    list(CR = df$gamma_selected, OLS = 0, CD = 1),
    function(gamma, method_name) {
      out_of_sample_risk_per_method(gamma) %>% mutate(method = method_name)
    })
  method_summary$method <- as.factor(method_summary$method)
  method_summary$in_sample_size <- as.numeric(
    as.character(method_summary$in_sample_size))
  method_summary$out_of_sample_shift <- factor(
    method_summary$out_of_sample_shift, levels = out_of_sample_shifts)

  ####################################################################
  # Plot performance of OLS, CR, and CD
  ####################################################################

  p1 <- ggplot(data = method_summary,
               mapping = aes(x = in_sample_size,
                             y = out_of_sample_risk_mean,
                             fill = method,
                             color = method)) +
    geom_line(alpha = 1) +
    geom_ribbon(mapping = aes(ymin = out_of_sample_risk_mean - out_of_sample_risk_sd,
                              ymax = out_of_sample_risk_mean + out_of_sample_risk_sd),
                alpha = 0.1,
                linewidth = 0) +
    facet_wrap(. ~ dim,
               ncol = length(levels(method_summary$dim)),
               scales = 'free_y',
               labeller = labeller(
                 dim = function(labels) {
                   paste("Dimension =", labels)
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
    scale_x_log10(
      labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    coord_cartesian(ylim = c(0, 0.5))

  save_plot(p1,
            paste0(path, "increasing_dim_comparison.pdf"),
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
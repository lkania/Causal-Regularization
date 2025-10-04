####################################################################
# Your R session's cwd should be the root of the project
####################################################################
source("./src/data.R")
source("./src/utils.R")
source("./src/measurements.R")
source("./src/selection.R")
source("./src/real.R")

library(digest)

# Note: This function must always return a data.frame
run_simulation <- function(trial,
                           target,
                           B,
                           out_of_sample_shift,
                           in_sample_size,
                           in_sample_shift,
                           intervention_mean,
                           in_sample_obs_mean,
                           folds_for_cv,
                           confounding,
                           shift_target) {

  out_of_sample_data <- gen(target = target,
                            mean = intervention_mean,
                            B = B,
                            n = out_of_sample_size,
                            shift_strength = out_of_sample_shift,
                            add_intercept = TRUE)

  in_sample_data <- data_gen(target = target,
                             B = B,
                             n = in_sample_size,
                             shift_strength = in_sample_shift,
                             obs_mean = in_sample_obs_mean,
                             intervention_mean = intervention_mean,
                             confounding = confounding,
                             shift_target = shift_target,
                             add_intercept = TRUE)

  results <- resample(
    data_train = in_sample_data,
    data_test = out_of_sample_data,
    compute_estimator_for_all_gamma = compute_estimator_for_all_gamma,
    folds_for_cv = folds_for_cv,
    metric_for_cv = metric_for_cv)

  results$out_of_sample_risk <- results$out_of_sample_risk / out_of_sample_shift
  results$in_sample_size <- in_sample_size
  results$in_sample_shift <- in_sample_shift
  results$in_sample_obs_mean <- in_sample_obs_mean
  results$intervention_mean <- intervention_mean
  results$out_of_sample_size <- out_of_sample_size
  results$out_of_sample_shift <- out_of_sample_shift
  results$B <- digest(B)
  results$trial <- trial

  return(results)

}
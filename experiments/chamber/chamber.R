####################################################################
# Your R session's working directory should be the root of the project
####################################################################
start_time <- Sys.time()

source("./src/parse.R")
folds_for_cv <- as.integer(get_arg_value("foldsForCV", 10))

source("./src/init_parallel.R")
source("./src/plotutils.R")
source("./src/real.R")
source("./src/measurements.R")
####################################################################
n_resamples <- 1000
compute_estimator_for_all_gamma <- (function(data) scr(data = data, length = 50))
metric_for_cv <- absolute_risk_difference
path <- './experiments/chamber/' # where to save the plot and find the data
####################################################################

variables <- c(
  "red",
  "green",
  "blue",
  "current",
  "pol_1",
  "pol_2",
  "ir_1",
  "ir_2",
  "ir_3",
  "vis_1",
  "vis_2",
  "vis_3",
  "l_11",
  "l_12",
  "l_21",
  "l_22",
  "l_31",
  "l_32"
)

covariates_to_use_1 <- c(
  "intercept",
  "red",
  "green",
  "blue"
)

covariates_to_use_2 <- c(
  "intercept",
  "red",
  "green",
  "blue",
  "current",
  "ir_1",
  "ir_3",
  "vis_1",
  "vis_2",
  "vis_3"
)

covariates_to_use_3 <- c(
  "intercept",
  "red",
  "green",
  "blue",
  "current",
  "pol_1",
  "pol_2",
  "ir_1",
  "ir_3",
  "vis_1",
  "vis_2",
  "vis_3",
  "l_11",
  "l_12",
  "l_21",
  "l_22",
  "l_31",
  "l_32"
)

get_dataset <- function(filenames) {
  data <- do.call(
    rbind,
    lapply(filenames, function(f) read.csv(paste0(path, 'data/', f), header = TRUE)))
  data <- data[, colnames(data) %in% variables]
  data$intercept <- 1
  data$y <- data$ir_2
  return(data)
}

data_ref <- get_dataset(filenames = c('uniform_reference.csv'))
data_int <- get_dataset(filenames = c(
  'uniform_red_mid.csv',
  'uniform_blue_mid.csv',
  'uniform_green_mid.csv'
))
data_out <- get_dataset(filenames = c(
  'uniform_red_strong.csv',
  'uniform_blue_strong.csv',
  'uniform_green_strong.csv'
))

#############################################
# Plot data
#############################################


# Plot ir_2 vs each color component
int_label <- 'In-sample interventional'
obs_label <- 'In-sample observational'
out_label <- 'Out of sample'
combined_data <- bind_rows(
  get_dataset(filenames = c('uniform_reference.csv')) %>% mutate(dataset = obs_label),
  get_dataset(filenames = c(
    'uniform_red_mid.csv')) %>% mutate(dataset = int_label),
  get_dataset(filenames = c(
    'uniform_red_strong.csv')) %>% mutate(dataset = out_label)
)

# Pivot to long format so we have color variable and its value
long_data <- combined_data %>%
  pivot_longer(cols = c(red, blue, green),
               names_to = "color",
               values_to = "color_value")

long_data$dataset <- factor(long_data$dataset,
                            levels = c(obs_label, int_label, out_label))

long_data <- long_data %>% filter(color == 'red')

p <- ggplot(long_data, aes(x = color_value,
                           y = ir_2,
                           linetype = dataset,
                           color = dataset)) +
  geom_point(alpha = 0.015) +
  geom_smooth(se = FALSE, alpha = 1) +
  labs(x = expression("Red Component Value " ~ "(" * R * ")"),
       y = expression("Target " ~ "(" * tilde(I)[2] * ")"),
       linetype = "Dataset",
       color = "Dataset") +
  theme(legend.position = "right",
        strip.text.x = element_text(size = 10),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12)) +
  guides(
    color = guide_legend(nrow = 3),
    linetype = guide_legend(nrow = 3)
  )


save_plot(p,
          paste0(path, "chamber_datasets.pdf"),
          dims = list(width = 6, height = 2))


#############################################
# Run the method
#############################################

set.seed(123)
source("./src/cr_l2.R")

run_resample_with_covariates <- function(covariates_to_use, label) {
  data_train <- list(
    Xe = as.matrix(data_int[, colnames(data_int) %in% covariates_to_use]),
    ye = data_int$y,
    Xo = as.matrix(data_ref[, colnames(data_ref) %in% covariates_to_use]),
    yo = data_ref$y
  )
  data_test <- list(
    X = as.matrix(data_out[, colnames(data_out) %in% covariates_to_use]),
    y = data_out$y
  )

  result <- resample(
    data_train = data_train,
    data_test = data_test,
    compute_estimator_for_all_gamma = compute_estimator_for_all_gamma,
    n_resamples = n_resamples,
    folds_for_cv = folds_for_cv,
    metric_for_cv = metric_for_cv
  )

  # Add label for identification
  result$covariate_set <- label
  return(result)
}

# Run and merge all results
result1 <- run_resample_with_covariates(covariates_to_use_1, "Only causal parameters")
result2 <- run_resample_with_covariates(covariates_to_use_2, "9 covariates")
result3 <- run_resample_with_covariates(covariates_to_use_3, "17 covariates")

# Merge all into a single data frame (assuming each `result` is a data frame)
df <- bind_rows(result1, result2, result3)
# Force a specific order of the factor for facet wrap
df$covariate_set <- factor(df$covariate_set, levels = c(
  "Only causal parameters",
  "9 covariates",
  "17 covariates"))

gamma_stats <- df %>%
  # there are many rows per idx, all of them have the same
  # gamma_selected and gamma_opt
  group_by(idx, covariate_set) %>%
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
  group_by(label, covariate_set) %>%
  summarise(
    gamma_mean = mean(gamma),
    gamma_sd = sd(gamma),
    gamma_lb = pmax(gamma_mean - gamma_sd, 0),
    gamma_ub = pmin(gamma_mean + gamma_sd, 1),
    .groups = "drop"
  )

p <- ggplot(df, aes(x = gammas, y = in_sample_risk_mean)) +
  geom_line(color = 'black') +
  geom_ribbon(mapping = aes(
    ymin = in_sample_risk_mean - in_sample_risk_sd,
    ymax = in_sample_risk_mean + in_sample_risk_sd),
              alpha = 0.1,
              linewidth = 0,
              fill = 'black') +
  geom_vline(data = gamma_stats %>% filter(label != "OPT"),
             mapping = aes(xintercept = gamma_mean, color = label),
             linetype = "dashed",
             size = 0.7) +
  geom_rect(data = gamma_stats %>% filter(label != "OPT"),
            mapping = aes(xmin = gamma_lb,
                          xmax = gamma_ub,
                          ymin = -Inf, ymax = Inf,
                          fill = label),
            alpha = 0.1,
            inherit.aes = FALSE) +
  facet_wrap(~covariate_set, scales = "free_y") +
  scale_x_continuous(
    breaks = seq(0, 1, 0.25),
    labels = expression(OLS, '0.25', '0.5', '0.75', CD)
  ) +
  scale_fill_discrete(name = "Method") +
  scale_color_discrete(name = "Method") +
  ylab('Cross-validated in-sample risk difference') +
  xlab(expression('Causal regularization path' ~ tilde(beta)[lambda])) +
  theme(legend.position = "top",
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  guides(color = guide_legend(nrow = 1))

save_plot(p,
          paste0(path, "chamber_in_sample_folds=", folds_for_cv, ".pdf"),
          dims = list(width = 10, height = 4))

p <- ggplot(df, aes(x = gammas, y = out_of_sample_risk_mean)) +
  geom_line(color = 'black') +
  geom_ribbon(mapping = aes(
    ymin = out_of_sample_risk_mean - out_of_sample_risk_sd,
    ymax = out_of_sample_risk_mean + out_of_sample_risk_sd),
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
  facet_wrap(~covariate_set, scales = "free_y") +
  scale_x_continuous(
    breaks = seq(0, 1, 0.25),
    labels = expression(OLS, '0.25', '0.5', '0.75', CD)
  ) +
  scale_fill_discrete(name = "Method") +
  scale_color_discrete(name = "Method") +
  ylab('Out-of-sample risk') +
  xlab(expression('Causal regularization path' ~ tilde(beta)[lambda])) +
  theme(legend.position = "top",
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  guides(color = guide_legend(nrow = 1))

save_plot(p,
          paste0(path, "chamber_out_of_sample_folds=", folds_for_cv, ".pdf"),
          dims = list(width = 10, height = 4))

####################################################################
# Print runtime
####################################################################
end_time <- Sys.time()
total_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
hours <- total_seconds %/% 3600
minutes <- (total_seconds %% 3600) %/% 60
seconds <- round(total_seconds %% 60)
cat(sprintf("Total runtime: %02d:%02d:%02d (hh:mm:ss)\n", hours, minutes, seconds))
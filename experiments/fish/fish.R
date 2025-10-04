####################################################################
# Your R session's working directory should be the root of the project
####################################################################
source("./src/init_parallel.R")
source("./src/utils.R")
source("./src/plotutils.R")
source("./src/real.R")
source("./src/measurements.R")
path <- './experiments/fish/'

####################################################################
n_resamples <- 1000
folds_for_cv <- 10
covariates_to_use <- c("intercept", "p")
compute_estimator_for_all_gamma <- (function(data) scr(data = data, length = 100))
metric_for_cv <- absolute_risk_difference
####################################################################

dataset <- loadRData(paste0(path, 'Fulton.Rdata'))
dataset$y <- dataset$q
dataset$intercept <- 1

mod1 <- dataset$Wed == 0
mod2 <- dataset$Stormy == 0
data_intervention_train_idx <- mod2 & mod1 # fair weather
data_obs_train_idx <- !mod2 & mod1
data_test_idx <- !mod1

data_intervention_train <- select(data = dataset,
                                  idx = data_intervention_train_idx,
                                  names = covariates_to_use)

data_obs_train <- select(data = dataset,
                         idx = data_obs_train_idx,
                         names = covariates_to_use)

data_train <- list(Xe = data_intervention_train$X,
                   ye = data_intervention_train$y,
                   Xo = data_obs_train$X,
                   yo = data_obs_train$y)

data_test <- select(data = dataset,
                    idx = data_test_idx,
                    names = covariates_to_use)

####################################################################
# visualize datasets
####################################################################

# visualize response
main_text_size <- 12
d_ <- dataset
d_$g <- "l"
d_[data_intervention_train_idx,]$g <- paste0(
  "In-Sample (Not Wednesday, Fair weather) (", sum(data_intervention_train_idx), " obs.)")
d_[data_obs_train_idx,]$g <- paste0(
  "In-Sample (Not Wednesday, Stormy weather) (", sum(data_obs_train_idx), " obs.)")
d_[data_test_idx,]$g <- paste0(
  "Out-of-Sample (Wednesday) (", sum(data_test_idx), " obs.)")
d_$p <- d_$p - mean(d_[data_obs_train_idx,]$p)
d_$y <- d_$y - mean(d_[data_obs_train_idx,]$y)

p0 <- ggplot(d_, aes(y = y, x = p, color = g)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  ylab('log(quantity)') +
  xlab('log(price)') +
  # scale_fill_discrete(name = "") +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(axis.text = element_text(size = main_text_size),
        legend.text = element_text(size = main_text_size))

p1 <- ggplot(d_, aes(y, group = g, fill = g)) +
  geom_density(alpha = 0.5) +
  ylab('Approximate density') +
  xlab('log(quantity)') +
  scale_fill_discrete(name = "") +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(axis.text = element_text(size = main_text_size),
        legend.text = element_text(size = main_text_size))

# visualize covariate
p2 <- ggplot(d_, aes(p, group = g, fill = g)) +
  geom_density(alpha = 0.5) +
  ylab('Approximate density') +
  xlab('log(price)') +
  scale_fill_discrete(name = "") +
  theme(legend.position = "top") +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(axis.text = element_text(size = main_text_size),
        legend.text = element_text(size = main_text_size))

pp2 <- ggarrange(p0, p2, p1, ncol = 3, nrow = 1, common.legend = TRUE, legend = "top")
save_plot(pp2, paste0(path, "fish_datasets.pdf"), dims = list(width = 15, height = 5))
####################################################################

set.seed(123)
source("./src/cr_l2.R")
invisible(resample(data_train = data_train,
                   data_test = data_test,
                   compute_estimator_for_all_gamma = compute_estimator_for_all_gamma,
                   n_resamples = n_resamples,
                   folds_for_cv = folds_for_cv,
                   metric_for_cv = metric_for_cv,
                   path_to_save_results = path,
                   prefix = 'fish'))
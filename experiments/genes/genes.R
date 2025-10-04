####################################################################
# Your R session's working directory should be the root of the project
####################################################################
source("./src/init_parallel.R")
source("./src/utils.R")
source("./src/plotutils.R")
source("./src/real.R")
path <- './experiments/genes/'

####################################################################
# Data description
# data$obs    data without interventions
# data$int    data with interventions
# data$intpos indices of the intervened genes
####################################################################
n_resamples <- 1000
n_covariates <- 10
folds_for_cv <- 10
compute_estimator_for_all_gamma <- (function(data) scr(data = data, length = 100))

data <- loadRData(paste0(path, 'Kemmeren.RData'))

intervened_idx <- as.vector(data$intpos)
intervened <- (1:6170 %in% intervened_idx)
not_intervened <- !intervened

writeLines(paste0("Number of observational samples ", dim(data$obs)[1]))
writeLines(paste0("Number of interventional samples ", dim(data$int)[1]))


# target selection
# Choose the most shifted non-intervened gene
observed_mean <- colMeans(data$obs) * as.numeric(not_intervened)
intervened_mean <- colMeans(data$int) * as.numeric(not_intervened)
diff <- abs(observed_mean - intervened_mean)
or <- order(diff, decreasing = TRUE)
target_idx <- or[1]
target_name <- names(diff[target_idx])[1]
writeLines(paste0("Gene ", target_idx, ' ', target_name,
                  " is the most shifted non-intervened gene. Shift ",
                  diff[target_idx]))

# covariate selection
# Choose the most shifted intervened genes
observed_mean <- colMeans(data$obs) * as.numeric(intervened)
intervened_mean <- colMeans(data$int) * as.numeric(intervened)
diff <- abs(observed_mean - intervened_mean)
or <- order(diff, decreasing = TRUE)

# most shifted intervened genes
covariates_idx <- or[c(1:n_covariates)]
covariates_names <- data$genenames[covariates_idx]
stopifnot(!(target_name %in% covariates_names))

yo <- as.matrix(data$obs[, target_idx])
stopifnot(!any(is.nan(yo)))

Xo <- as.matrix(data$obs[, covariates_idx])
stopifnot(!any(is.nan(Xo)))

ye <- as.matrix(data$int[, target_idx])
stopifnot(!any(is.nan(ye)))

Xe <- as.matrix(data$int[, covariates_idx])
stopifnot(!any(is.nan(Xe)))

half <- floor(length(ye) / 2)
idx <- ((1:length(ye)) <= half)

data_train <- list(Xe = as.matrix(Xe[idx,]),
                   ye = as.vector(ye[idx]),
                   Xo = as.matrix(Xo),
                   yo = as.vector(yo))


data_train$Xe <- cbind(data_train$Xe, intercept = rep(1, length(data_train$ye)))
data_train$Xo <- cbind(data_train$Xo, intercept = rep(1, length(data_train$yo)))

data_test <- list(X = as.matrix(Xe[!idx,]),
                  y = as.vector(ye[!idx]))
#add intercept
data_test$X <- cbind(data_test$X, rep(1, intercept = length(data_test$y)))


set.seed(123)
source("./src/cr_l2.R")
invisible(resample(data_train = data_train,
                   data_test = data_test,
                   compute_estimator_for_all_gamma = compute_estimator_for_all_gamma,
                   n_resamples = n_resamples,
                   folds_for_cv = folds_for_cv,
                   metric_for_cv = absolute_risk_difference,
                   path_to_save_results = path,
                   prefix = 'genes'))
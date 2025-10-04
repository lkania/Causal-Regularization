library(future.apply)
# Detect the number of available logical CPU cores
cores <- availableCores()

# Determine number of cores to be used
n_jobs <- switch(as.character(cores),
                 "32" = 20,
                 "12" = 10,
                 "4" = 2,
                 "1" = 1,
                 as.integer(cores / 2)
)

# Set the plan with the determined number of workers
writeLines(sprintf("Using %d of %d available cores.\n", n_jobs, cores))

# Note: if you use Windows or RStudio, you should use:
#   plan(multisession, workers = n_jobs)
# instead of the command below
plan(multicore, workers = n_jobs)
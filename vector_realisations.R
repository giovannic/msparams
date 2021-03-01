library(parallel)

args = commandArgs(trailingOnly=TRUE)
node <- args[1]
in_dir <- args[2]
out_dir <- args[3]
lib <- args[4]
batch_size <- 50

print(paste0('beginning node ', node))

params <- read.csv(file.path(in_dir, paste0('p_space_vector_', node, '.csv')))

# setup cluster...
cl <- makeCluster(detectCores())
clusterExport(cl, "params")
clusterExport(cl, "lib")

# set up processing function
clusterEvalQ(cl, {
  library(malariasimulation, lib=lib)

  population <- 1000

  process_row <- function(row, p_row) {
    year <- 365
    month <- 30
    sim_length <- 10 * year

    simparams <- get_parameters(
      list(
        human_population = population,
        average_age = row$average_age,
        model_seasonality = TRUE,
        g0 = row$seasonal_a0,
        g = as.numeric(row[c('seasonal_a1', 'seasonal_a2', 'seasonal_a3')]),
        h = as.numeric(row[c('seasonal_b1', 'seasonal_b2', 'seasonal_b3')]),
        Q0 = row$Q0
      )
    )

    simparams <- parameterise_total_M(simparams, population * row$total_M)

    output <- run_simulation(
      sim_length,
      simparams
    )[c('timestep', 'pv_730_3650')]

    # take last 5 years
    output <- tail(output, 5 * year)

    output[['p_row']] <- p_row
    output
  }
})

batches <- split(
  seq(nrow(params)),
  (seq(nrow(params))-1) %/% batch_size
)

for (batch_i in seq_along(batches)) {
  start_time <- Sys.time()
  print(paste0('node ', node, ' batch ', batch_i, ' starting'))
  # do the work
  results <- do.call(
    'rbind',
    parLapply(
      cl,
      X = batches[[batch_i]],
      f = function (i) process_row(params[i,], i)
    )
  )

  write.csv(
    results,
    file.path(out_dir, paste0('realisation_', node, '_batch_', batch_i, '.csv')),
    row.names = FALSE
  )
  print(paste0('node ', node, ' batch ', batch_i, ' completed'))
  print(Sys.time())
  print(Sys.time() - start_time)
}


stopCluster(cl)


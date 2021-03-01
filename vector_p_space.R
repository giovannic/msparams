library(reshape2)

out_dir <- './msparams/vector_p_space_v1/'

n_realisations <- 5e+5

set.seed(42)
print('runtime (days)')
runtime <- 4 #mins
cores <- 8
nodes <- 32
print(n_realisations * runtime / (60 * 24 * cores * nodes))

history <- readRDS('./msparams/GTS2020_sites_fitted.RDS')
seasonality <- do.call(rbind, history$seasonality)

vector_params <- data.frame(
  average_age = runif(n_realisations, 20 * 365, 40 * 365),
  total_M = runif(n_realisations, 0, 100),
  Q0 = runif(n_realisations, 0, 1)
)

vector_params <- cbind(
  vector_params,
  seasonality[sample(nrow(seasonality), n_realisations, replace = TRUE), ]
)

node_map <- rep(1:nodes, each=nodes, length.out=n_realisations)

write_out_splits <- function(splits, name) {
  for (i in seq(nodes)) {
    write.csv(
      splits[[i]],
      file.path(out_dir, paste0('p_space_', name, '_', i, '.csv')),
      row.names = FALSE
    )
  }
}

write_out_splits(split(vector_params, node_map), 'vector')

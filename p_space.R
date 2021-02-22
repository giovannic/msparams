library(reshape2)

out_dir <- './msparams/p_space_v3/'

n_realisations <- 5e+5

print('runtime (days)')
runtime <- 4 #mins
cores <- 8
nodes <- 32
print(n_realisations * runtime / (60 * 24 * cores * nodes))

history <- readRDS('./msparams/GTS2020_sites_fitted.RDS')
seasonality <- do.call(rbind, history$seasonality)

vector_profile <- function (r) {
  props <- r[c('prop_gamb_ss', 'prop_fun', 'prop_arab')]
  c(
    Q0 = weighted.mean(r[c('gamb_ss_Q0', 'fun_Q0', 'arab_Q0')], props),
    rs = weighted.mean(r[c('irs_dif_gamb_ss', 'irs_dif_fun', 'irs_dif_arab')], props), # is this right?
    phi_bednets = weighted.mean(r[c('gamb_ss_Q_bed', 'fun_Q_bed', 'arab_Q_bed')], props),
    phi_spraying = weighted.mean(r[c('gamb_ss_Q_in', 'fun_Q_in', 'arab_Q_in')], props)
  )
}

vector_profiles <- do.call(rbind, lapply(history$vectors, vector_profile))

create_interventions <- function (row) {
  cbind(row, history$interventions[[row]])
}

cast_intervention_history <- function(interventions, columns) {
  do.call(cbind, lapply(columns, function(colname) {
    wide <- dcast(interventions, row ~ year, value.var=colname, drop=FALSE)
    wide <- wide[, names(wide) != 'row']
    colnames(wide) <- paste0(colname, '_', colnames(wide))
    wide
  }))
}

interventions <- do.call(rbind, lapply(seq(nrow(history)), create_interventions))

synthetic_params <- data.frame(
  average_age = runif(n_realisations, 20 * 365, 40 * 365),
  ft = runif(n_realisations, 0, 1),
  prop_act = runif(n_realisations, 0, 1),
  vaccine_rounds_per_year = sample(seq(from = 0, to = 1), n_realisations, replace = TRUE),
  mda_rounds_per_year = sample(seq(from = 0, to = 2), n_realisations, replace = TRUE),
  smc_rounds_per_year = sample(seq(from = 0, to = 2), n_realisations, replace = TRUE),
  net_rounds_per_year = sample(seq(from = 0, to = 2), n_realisations, replace = TRUE),
  spraying_rounds_per_year = sample(seq(from = 0, to = 2), n_realisations, replace = TRUE), # one or two (perennial) rounds per year before the peak
  vaccine_smc_correlation = runif(n_realisations, -1, 1),
  vaccine_net_correlation = runif(n_realisations, -1, 1),
  spraying_net_correlation = runif(n_realisations, -1, 1)
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

write_out <- function(data, name) {
  write.csv(
    data,
    file.path(out_dir, paste0('p_space_', name, '.csv')),
    row.names = FALSE
  )
}

write_out_splits(split(synthetic_params, node_map), 'synthetic')
write_out_splits(split(sample.int(dim(history)[[1]], n_realisations, replace=TRUE), node_map), 'location')
write_out(seasonality, 'season')
write_out(as.data.frame(vector_profiles), 'vectors')

write_out(
  history[c('Continent', 'NAME_0', 'NAME_1', 'ur', 'total_M')],
  'basics'
)

i_names <- c('tx', 'prop_act', 'llin', 'irs', 'irs_rounds', 'smc_rounds', 'rtss')
for (i_name in i_names) {
  write_out(
    cast_intervention_history(interventions, i_name),
    i_name
  )
}

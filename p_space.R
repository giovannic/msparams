args = commandArgs(trailingOnly=TRUE)
out_dir <- args[1]

p_space = list(
  human_population = seq(from = 5000, to = 9000, by = 500),
  average_age = seq(from = 20 * 365, to = 40 * 365, by = 5 * 365),
  seasonality_profile = c('bimodal', 'high', 'perennial', 'normal'),
  eir = c(.1, .5, 1, 10, 50, 75, 100),
  species_12 = seq(from = 0, to = 1, by = .2),
  species_23 = seq(from = 0, to = 1, by = .2),
  ft = seq(from = 0, to = 1, by = .2),
  drug_balance = seq(from = 0, to = 1, by = .2),
  vaccine_start = seq(from = 0, to = 4) * 365,
  vaccine_duration = seq(5) * 365,
  vaccine_frequency = seq(from = .5, to = 2, by = .5) * 365,
  vaccine_coverage = seq(from = 0, to = 1, by = .2),
  booster_delay = seq(from = .5, to = 2, by = .5) * 365,
  booster_coverage = seq(from = 0, to = 1, by = .2),
  mda_start = seq(from = 0, to = 4) * 365,
  mda_duration = seq(5) * 365,
  mda_frequency = seq(from = .5, to = 2, by = .5) * 365,
  mda_coverage = seq(from = 0, to = 1, by = .2),
  smc_start = seq(from = 0, to = 4) * 365,
  smc_duration = seq(5) * 365,
  smc_frequency = seq(from = .5, to = 2, by = .5) * 365,
  smc_coverage = seq(from = 0, to = 1, by = .2),
  net_start = seq(from = 0, to = 4) * 365,
  net_duration = seq(5) * 365,
  net_coverage = seq(from = 0, to = 1, by = .2),
  net_retention = seq(5) * 365,
  spraying_start = seq(from = 0, to = 4) * 365,
  spraying_duration = seq(5) * 365,
  spraying_coverage = seq(from = 0, to = 1, by = .2),
  rtss_smc_correlated = c(1, 0, -1),
  rtss_net_correlation = c(1, 0, -1)
)

p_space = list(
  human_population = seq(from = 5000, to = 9000, by = 500),
  average_age = seq(from = 20 * 365, to = 40 * 365, by = 10 * 365),
  seasonality_profile = c('bimodal', 'high', 'perennial', 'normal'),
  eir = c(.1, .5, 1, 10, 50, 75, 100),
  ft = seq(from = 0, to = 1, by = .2),
  vaccine_frequency = seq(from = .5, to = 1, by = .5) * 365,
  vaccine_coverage = c(0, .6, .8),
  smc_frequency = seq(from = .5, to = 1, by = .5) * 365,
  smc_coverage = c(0, .6, .8),
  net_frequency = seq(from = .5, to = 1, by = .5) * 365,
  net_coverage = c(0, .6, .8)
)

# Remove all interventions (except vaccines)
# Remove boosters
# Remove species
# Simple seasonality
# Remove correlation
# No repetitions

print('dims')
print(length(p_space))
print('dim counts')
counts <- vapply(p_space, function(s) length(s), numeric(1))
names(counts) <- names(p_space)
print(counts)
print('realisations')
print(prod(counts))
print('days')
runtime <- 1 #mins
cores <- 8
nodes <- 32
print(prod(counts) * runtime / (60 * 24 * cores * nodes))

params <- expand.grid(p_space)

nr <- nrow(params)
outs <- split(params, rep(1:nodes, each=nodes, length.out=nr))

for (i in seq_along(outs)) {
  write.csv(
    outs[[i]],
    file.path(out_dir, paste0('p_space_', i, '.csv')),
    row.names = FALSE
  )
}
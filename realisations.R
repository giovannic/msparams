library(parallel)

args = commandArgs(trailingOnly=TRUE)
node <- args[1]
in_dir <- args[2]
out_dir <- args[3]
lib <- args[4]
batch_size <- 50

print(paste0('beginning node ', node))

interventions <- cbind(
  read.csv(file.path(in_dir, paste0('p_space_synthetic_', node, '.csv'))),
  read.csv(file.path(in_dir, paste0('p_space_location_', node, '.csv')))
)

locations <- cbind(
  read.csv(file.path(in_dir, paste0('p_space_season.csv'))),
  read.csv(file.path(in_dir, paste0('p_space_vectors.csv'))),
  read.csv(file.path(in_dir, paste0('p_space_basics.csv')))
)

params <- merge(interventions, locations, by.x = 'x', by.y = 0)

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
    sim_length <- 5 * year
    
    simparams <- get_parameters(
      list(
        human_population = population,
        variety_proportions = 1,
        average_age = row$average_age,
        total_M = ceiling(population * row$total_M),
        mosquito_limit = ceiling(population * row$total_M * 10),
        model_seasonality = TRUE,
        g0 = row$seasonal_a0,
        g = as.numeric(row[c('seasonal_a1', 'seasonal_a2', 'seasonal_a3')]),
        h = as.numeric(row[c('seasonal_b1', 'seasonal_b2', 'seasonal_b3')]),
        Q0 = row$Q0,
        rs = row$rs,
        phi_bednets = row$phi_bednets,
        phi_indoors = row$phi_spraying,
        endophily = .813,
        blood_meal_rates = 1/3,
        rn = .56,
        rnm = .24,
        dn0 = .533
      )
    )
    
    #TODO: human equilibrium from Total_M
    simparams <- set_drugs(simparams, list(AL_params, DHC_PQP_params, SP_AQ_params))
    simparams <- set_clinical_treatment(
      simparams,
      row$ft,
      c(1, 2),
      c(1 - row$prop_act, row$prop_act)
    )
    
    peak <- peak_season_offset(simparams)
    
    if (row$vaccine_rounds_per_year > 0) {
      simparams <- set_rtss(
        simparams,
        start = 0,
        end = sim_length,
        frequency = floor(year / row$vaccine_rounds_per_year),
        coverage = .9,
        min_ages = 5 * month,
        max_ages = 17 * month,
        boosters = 18 * month,
        booster_coverage = .7
      )
    }

    if (row$net_rounds_per_year > 0) {
      simparams <- set_bednets(
        simparams,
        timesteps = seq(10) * floor(year / row$net_rounds_per_year) + peak - 3 * month,
        coverages = rep(.9, 10),
        retention = 5 * year
      )
    }
    
    if (row$mda_rounds_per_year > 0) {
      simparams <- set_mda(
        simparams,
        drug = 3,
        start = 0,
        end = sim_length,
        frequency = floor(year / row$mda_rounds_per_year),
        min_age = 0,
        max_age = 200 * year,
        .9
      )
    }
    
    if (row$smc_rounds_per_year > 0) {
      simparams <- set_smc(
        simparams,
        drug = 3,
        start = 0 + peak - 3 * month,
        end = sim_length,
        frequency = floor(year / row$smc_rounds_per_year),
        min_age = 2 * year - 1,
        max_age = 11 * year,
        .9
      )
    }
    
    if (row$spraying_rounds_per_year > 0) {
      simparams <- set_spraying(
        simparams,
        timesteps = seq(10) * floor(year / row$spraying_rounds_per_year) + peak - 3 * month,
        coverages = rep(.9, 10)
      )
    }
    
    correlations = get_correlation_parameters(simparams)
    if (row$vaccine_rounds_per_year > 0 && row$smc_rounds_per_year > 0) {
      correlations$inter_intervention_rho('rtss', 'smc', row$vaccine_smc_correlation)
    }
    if (row$vaccine_rounds_per_year > 0 && row$net_rounds_per_year > 0) {
      correlations$inter_intervention_rho('rtss', 'bednets', row$vaccine_net_correlation)
    }
    if (row$net_rounds_per_year > 0 && row$spraying_rounds_per_year > 0) {
      correlations$inter_intervention_rho('bednets', 'spraying', row$spraying_net_correlation)
    }
    
    set.seed(42)
    output <- run_simulation(
      sim_length,
      simparams,
      correlations = correlations
    )[c('timestep', 'pv_730_3650')]
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


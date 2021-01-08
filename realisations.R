library(parallel)

args = commandArgs(trailingOnly=TRUE)
node <- args[1]
in_dir <- args[2]
out_dir <- args[3]

p_space <- file.path(in_dir, paste0('p_space_', node, '.csv'))
params <- head(read.csv(p_space), 4)

# setup cluster...
numCores <- detectCores()
cl <- makeCluster(numCores)
clusterExport(cl, "params")

# set up processing function
clusterEvalQ(cl, {
  library(malariasimulation)
  library(malariaEquilibrium)
  
  seasonality <- list(
    bimodal = list(
      g0 = 0.28605,
      g = c(0.20636, -0.0740318, -0.0009293),
      h = c(0.173743, -0.0730962, -0.116019)
    ),
    high = list(
      g0 = 0.284596,
      g = c(-0.317878, -0.0017527, 0.116455),
      h	= c(-0.331361, 0.293128, -0.0617547)
    ),
    perennial = list(
      g0 = 0.285277,
      g = c(-0.0248801, -0.0529426, -0.016891),
      h = c(-0.0216681, -0.0242904, -0.0073646)
    ),
    normal = list(
      g0 = 0.285505,
      g = c(-0.325352, -0.0109352, 0.0779865),
      h =	c(-0.132815, 0.104675, -0.013919)
    )
  )
  
  process_row <- function(row, p_row) {
    year <- 365
    month <- 30
    sim_length <- 5 * year
    jamie_params <- load_parameter_set()
    
    simparams <- get_parameters(c(
      translate_jamie(remove_unused_jamie(jamie_params)),
      list(
        human_population = row$human_population,
        variety_proportions = 1,
        model_seasonality = TRUE,
        average_age = row$average_age
      ),
      seasonality[[row$seasonality_profile]]
    ))
    eq <- human_equilibrium(EIR = row$eir, ft = row$ft, p = jamie_params, age = 0:99)
    simparams <- parameterise_human_equilibrium(simparams, eq)
    simparams <- parameterise_mosquito_equilibrium(simparams, EIR=row$eir)
    simparams <- set_drugs(simparams, list(AL_params, DHC_PQP_params))
    simparams <- set_clinical_treatment(simparams, row$ft, c(1, 2), c(.5, .5))
    simparams <- set_rtss(
      simparams,
      start = 0,
      end = sim_length,
      frequency = row$vaccine_frequency,
      coverage = row$vaccine_coverage,
      min_ages = 5 * month,
      max_ages = 17 * month,
      boosters = 18 * month,
      booster_coverage = .7
    )
    simparams <- set_bednets(
      simparams,
      timesteps = seq(10) * row$net_frequency * year,
      coverages = rep(row$net_coverage, 10),
      retention = 5 * year
    )
    
    output <- run_simulation(sim_length, simparams)[c('timestep', 'pv_730_3650')]
    output[['p_row']] <- p_row
    output
  }
})

# do the work
results <- do.call(
  'rbind',
  parLapply(
    cl,
    X = seq(nrow(params)),
    f = function (i) process_row(params[i,], i)
  )
)

stopCluster(cl)

write.csv(
  results,
  file.path(out_dir, paste0('realisation_', node, '.csv')),
  row.names = FALSE
)


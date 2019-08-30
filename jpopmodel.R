# R script to implement a PVA model for Brazilian Jaguar in the Atlantic Forest region.
# Author: Ascelin Gordon ascelin.gordon@rmit.edu.au
# Commencement Date: 02 Feb 2016

# To run:
#   source( 'jpopmodel.R')

rm(list=ls(all=TRUE))

#load required libraries
library(triangle)
library(actuar)

# set seed if want to be able to reproduce results
set.seed(44)

#source files
source('jpopmodel.functions.R')
source('jpopmodel.main.R')


        # ------------------------------------------------
        # Define parameters of the model
        # ------------------------------------------------

time.steps <- 20
stochastic.realizations <- 5
expert.realizations <- 3

use.fast.simulated.data <- FALSE  # if true use a small set of JCUs for fast testing

output.filename <- 'jpopmodel_data_disp.Rdata'


        # ------------------------------------------------
        # Run the model
        # ------------------------------------------------

par(mfrow=c(2,3))

if(!use.fast.simulated.data) {

  # in this case sampling from the elicited distributions

  # read in the elicited parameters
  ParamElic <- read.csv('input_data/elicited_parameters.csv')
  Distances <- read.csv('input_data/distmatrix.csv',header=F)
  JCUData <- read.csv('input_data/jcu_hunting_start_pop_unocc.csv')

  # Note if use this file, only has 18 JCUs (no Unocc JCUs, so code runs much faster)
  #JCUData <- read.csv('input_data/jcu_hunting_start_pop.csv')

  # sample parameter values based on elicited values
  Params <- get.parameters(Elicitation = ParamElic, JCUs = JCUData, Distances = Distances, Reps = expert.realizations)

  # run population model
  model.output <- apply(X = Params$ENSEMBLE.LIST, MARGIN = 1, FUN = run.pop.model.apply, Params.List = Params, Years = time.steps, Reps = stochastic.realizations)

  # convert the list of arrays that is an output from apply to a single array by binding them by row
  model.output <-  do.call(rbind, model.output)

} else {

  # generate some data for doing fast tests, 
  # currently set to havge just 3 JCUs
  source('gen.test.data.R')

  # create initial output data for fist run. Assume just one realization per expert
  cat('\nExpert 1')
  model.output <- run.jpop.model(expert.ID=1, expert.realization=1, stochastic.realizations, initial.population, jcu.attributes, disp.mort.matrix, time.steps)

  # dummy loop through experts
  for( expert in 2:3 ) {
    cat('\nExpert', expert)
    x <- run.jpop.model(expert, expert.realization=1, stochastic.realizations, initial.population, jcu.attributes, disp.mort.matrix, time.steps  )

    # append the run the current outputs
    model.output <- rbind( model.output, x)
  }

}

        # ------------------------------------------------
        # Save the outputs for later analysis
        # ------------------------------------------------


# Save the outputs to file. This is used in the analyse.results.R
# script. Note this currently gets overwritten each time the script
# runs.

saveRDS(model.output, file=output.filename)
cat('\nWrote ', output.filename, '\n')

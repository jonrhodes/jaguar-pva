# R script to implement a PVA model for Brazilian Jaguar in the Atlantic Forest region.
# Author: Ascelin Gordon ascelin.gordon@rmit.edu.au
# Commencement Date: 02 Feb 2016

# To run:
#   source( 'jpopmodel.R')

rm(list=ls(all=TRUE))

#load required libraries
library(triangle)

# set seed if want to be able to reproduce results
set.seed(44)

#source files
source('jpopmodel.functions.R')
source('jpopmodel.main.R')

        # ------------------------------------------------
        # Define parameters of the model
        # ------------------------------------------------

time.steps <- 50
stochastic.realizations <- 2
expert.realizations <- 2
output.filename <- 'jpopmodel_data_disp.Rdata'

#elicited parameters
ParamElic <- read.csv('input_data/elicited_parameters.csv')
Distances <- read.csv('input_data/distmatrix.csv',header=F)
JCUData <- read.csv('input_data/jcu_hunting_start_pop.csv')

#get random parameter values based on elicited values
Params <- get.parameters(Elicitation = ParamElic,JCUs = JCUData,Distances = Distances,Reps = expert.realizations)

#run population model
Output <- apply(X = Params$ENSEMBLE.LIST, MARGIN = 1, FUN = run.pop.model.apply, Params.List = Params, Years = time.steps, Reps = stochastic.realizations)

# generate some data for doing fast tests
source('gen.test.data.R')


        # ------------------------------------------------
        # Run the model
        # ------------------------------------------------

# create initial output data for fist run
cat('\nExpert 1')
model.output <- run.jpop.model(expert.ID=1, expert.realization, stochastic.realizations, initial.population, jcu.attributes, disp.mort.matrix, time.steps)

# dummy loop through experts
for( expert in 2:3 ) {
 cat('\nExpert', expert)
  #do a subsequent run to a append
  x <- run.jpop.model(expert, expert.realization, stochastic.realizations, initial.population, jcu.attributes, disp.mort.matrix, time.steps  )

  #append the run the current outputs 
  model.output <- rbind( model.output, x)

}


        # ------------------------------------------------
        # Save the outputs for later analysis
        # ------------------------------------------------


# Save the outputs to file. This is used in the analyse.results.R
# script. Note this currently gets overwritten each time the script
# runs.

saveRDS(model.output, file=output.filename)
cat('\nWrote ', output.filename, '\n')


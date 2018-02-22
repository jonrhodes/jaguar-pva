# R script to implement a PVA model for Brazilian Jaguar in the Atlandtic Forest region.
# Author: Ascelin Gordon ascelin.gordon@rmit.edu.au
# Commencement Date: 02 Feb 2016


# To run:
#   source( 'jpopmodel.R')

rm(list=ls(all=TRUE))

# set seed if want to be able to reproduce results
set.seed(2)

source('jpopmodel.functions.R')
source('jpopmodel.main.R')


        # ------------------------------------------------
        # Define parameters of the model
        # ------------------------------------------------


time.steps <- 40
output.filename <- 'jpopmodel_data_disp.Rdata'
stochastic.realizations <- 50

# Jonathan's code to provide, setting to dummy value for now.
expert.ID = -999
expert.realization = -999


# generate some data for doing fast tests
source('gen.test.data.R')

model.output <- run.jpop.model(expert.ID, expert.realization, stochastic.realizations, initial.population, jcu.attributes, disp.mort.matrix, time.steps  )

# Save the outputs to file. This is used in the analyse.results.R
# script. Note this currently gets overwritten each time the script
# runs.

saveRDS(model.output, file=output.filename)
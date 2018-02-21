# R script to implement a PVA model for Brazilian Jaguar in the Atlandtic Forest region.
# Author: Ascelin Gordon ascelin.gordon@rmit.edu.au
# Commencement Date: 02 Feb 2016


# To run:
#   source( 'jpopmodel.R')

rm(list=ls(all=TRUE))

set.seed(2)

source('jpopmodel.functions.R')
source('jpopmodel.main.R')

        # ------------------------------------------------
        # Define the name parameters of the model
        # ------------------------------------------------

DEBUG.LEVEL <- 0  # 0-none; 1-terse, 2-verbose
# num.time.steps <- 80
# num.reps <- 20
num.time.steps <- 40
num.reps <- 20
OPT.INCLUDE.DISPERSAL <- TRUE
output.filename <- 'jpopmodel_data_disp.Rdata'

num.life.stages <- 3

# Option of whether to use a truncated Poisson distribution for determing the
# number offspring (TRUE) or whether to assume there are only 1 or 2 offspring
# (FALSE). The truncated Poisson means that occasionally there may be 5 or 6
# offspring.
OPT.USE.TRUNCATED.POISSON.FOR.REPRODUCTION <- TRUE 


source('gen.test.data.R')

run.jpop.model(initial.population, jcu.attributes, disp.mort.matrix)
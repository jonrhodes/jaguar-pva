

# R script to implement a PVA model for Brazilian Jaguar... 
# Author: Ascelin Gordon ascelin.gordon@rmit.edu.au
# Commencement Date: 02 Feb 2016


# To run:
#   source( 'jpopmodel.R')

rm(list=ls(all=TRUE))


source('jpopmodel.functions.R')



DEBUG <- FALSE

# initialize variables
num.time.steps <- 40
num.reps <- 5


        # ------------------------------------------------
        # These are to be automatically extracted by the R
        # preprocessing scripts that Hawthorne will is writing
        # ------------------------------------------------


num.jcus <- 4
num.life.stages <- 3
jcu.att <- data.frame( cc=c(10,20,7,50),
                      #death.rate=c(0.5,0.01,0.002,0.4) #assuming same for each stage for now so one no per JCU
                      death.rate=rep(0.2,num.jcus),  # assuming same for each stage for now so one no per JCU
                      birth.rate=rep(0.2, num.jcus) # assuming only last stage gives birth
                      )


# There are 4 age classes (0-1 year olds, 1-2 year olds, 2-3 year olds, and 3+ year olds).

# make an initial population (random from now) 
initial.pop <- matrix( ncol=num.jcus, nrow=num.life.stages, sample(1:10, size=num.jcus*3, replace=TRUE ) )

colnames(initial.pop) <- paste( 'jcu', 1:num.jcus, sep='') # set column names
rownames(initial.pop) <- paste( 'stage', 1:num.life.stages, sep='') # set column names



        # ------------------------------------------------
        # Do validity checks
        # ------------------------------------------------


if( dim(jcu.att)[1] != num.jcus ) stop( '\n\nERROR JCU attributes not conistent with number of JCUs')

        # ------------------------------------------------
        # Make some objects to store results as model runs
        #
        # ------------------------------------------------

# a data frame to store the full state of the system
all.outputs <- data.frame ( "rep"=NA, "time"=NA, "jcu"=NA,
                           "stage1"=NA, "stage2"=NA, "stage3"=NA ) [numeric(0), ]

# a datafram to just store the total population size to make a quick plot at the end
total.pop <- matrix( ncol = num.time.steps, nrow = num.reps )

    
        # ------------------------------------------------
        # Do the model realizations 
        # ------------------------------------------------



for( rep in 1:num.reps) {
    
    if(DEBUG){cat( "\n\n--------------------------------"); cat( '\n Rep = ', rep, '\n' )}

    current.pop <- initial.pop
    
    
    for( time in 1:num.time.steps) {


        
        if( time > 1 ) {            
            # age the population
            current.pop <- age.population(current.pop, num.life.stages)
            # do reproduction (assumes for now only last stage reporduces)
            current.pop <- reproduce(current.pop, num.life.stages, jcu.att$birth.rate)

        }
        
        # loop through each JCU and apply the JCU specific mortality to it
        for( jcu in 1:num.jcus){

            if( time > 1 ) { # don't apply mortality in the first time step
                # Apply mortality to the populations in each JCU
                #y = function(x) x+1 # pop of each stage increases by 1
                y = function(x) rbinom(n=1, size=x, prob=(1-jcu.att$death.rate[jcu]))
                current.pop[,jcu] <- sapply( current.pop[,jcu], y)
            }
                
            # build a vector for the current info of the system
            current.info <- c(rep, time, jcu, current.pop[,jcu] )
                
            # save all the current info
            if( dim(all.outputs)[1] == 0 )
                # if first data entry, replace the first line
                all.outputs[1,] <- current.info
            else
                # otherwise append it
                all.outputs <- rbind( all.outputs, current.info )
        }


        total.pop[rep, time] <- sum(current.pop)

        #browser()
        
        if(DEBUG){ cat( '\n Time step = ', time, '\n' ); show(current.pop) }
        
    }
    
    if(DEBUG) cat("\n")

  
}

# plot the total pop trajectory of each realisationx and also the mean
# trajectory
matplot ( t(total.pop), type = 'l', main= 'Total population', ylab='pop size'  )
mean.traj <- apply(total.pop,2,mean)
lines(mean.traj, lwd=2)


save(all.outputs, file='jpopmodel_data.Rdata')

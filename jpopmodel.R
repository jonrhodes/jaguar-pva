

# R script to implement a PVA model for Brazilian Jaguar... 
# Author: Ascelin Gordon ascelin.gordon@rmit.edu.au
# Commencement Date: 02 Feb 2016


# To run:
#   source( 'jpopmodel.R')

rm(list=ls(all=TRUE))


source('jpopmodel.functions.R')


        # ------------------------------------------------
        # Define the name parameters of the model
        # ------------------------------------------------

DEBUG <- FALSE
num.time.steps <- 50
num.reps <- 10


        # ------------------------------------------------
        # It's planned that These parameters will be eventually be
        # extracted from the GIS and other data data we will be
        # given. But just specifying them here for now
        # ------------------------------------------------


num.jcus <- 4
num.life.stages <- 3

jcu.att <- data.frame( cc=c(10,20,7,50),
                      #mortality=c(0.5,0.1,0.2,0.4), #assuming same for each stage for now so one no per JCU
                      mortality=rep(0.2,num.jcus),
                      
                      # Assuming only last stage gives birth, this is the prob of giving birth to 1 offspring
                      birth.rate=rep(0.2, num.jcus)
                      
                      )

# Make an initial population (random from now) 
initial.pop <- matrix( ncol=num.jcus, nrow=num.life.stages, sample(1:10, size=num.jcus*3, replace=TRUE ) )
colnames(initial.pop) <- paste( 'jcu', 1:num.jcus, sep='') # set column names
rownames(initial.pop) <- paste( 'stage', 1:num.life.stages, sep='') # set the row names


        # ------------------------------------------------
        # Validity checks
        # ------------------------------------------------

if( dim(jcu.att)[1] != num.jcus ) stop( '\n\nERROR JCU attributes not conistent with number of JCUs')


        # ------------------------------------------------
        # Make some objects to store results as the model runs
        # ------------------------------------------------

# A data.frame to store the full state of the system
all.outputs <- data.frame ( "rep"=NA, "time"=NA, "jcu"=NA,
                           "stage1"=NA, "stage2"=NA, "stage3"=NA ) [numeric(0), ]

# A datafram to just store the total population size to make a quick plot at the end
total.pop <- matrix( ncol = num.time.steps, nrow = num.reps )

    
        # ------------------------------------------------
        # Run the model
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
        
        # loop through each JCU and apply the JCU specific mortality to each life stage
        for( jcu in 1:num.jcus){
            
            if( time > 1 ) { # don't apply mortality in the first time step
                # Apply mortality to the populations in each JCU
                current.pop[,jcu] <- apply.mortality(current.pop[,jcu], jcu.att$mortality[jcu])                
            }

            
            # ------------------------------------------------
            # Save the current info 
            # ------------------------------------------------
                
            # Build a vector for the current info of the system
            current.info <- c(rep, time, jcu, current.pop[,jcu] )
                
            if( dim(all.outputs)[1] == 0 )
                # if first data entry, replace the first line
                all.outputs[1,] <- current.info
            else
                # otherwise append it
                all.outputs <- rbind( all.outputs, current.info )
        }
        
        # Also save the total population size to make a quick plot 
        total.pop[rep, time] <- sum(current.pop)

        #browser()
        
        if(DEBUG){ cat( '\n Time step = ', time, '\n' ); show(current.pop) }
        
    }
    
    if(DEBUG) cat("\n")

  
}


        # ------------------------------------------------
        # Make some simple plots. For more detailed analysis run the
        # analyse.results.R script
        # ------------------------------------------------


# Plot the total pop trajectory of each realisationx and also the mean
# trajectory
matplot ( t(total.pop), type = 'l', main= 'Total population', ylab='pop size'  )
mean.traj <- apply(total.pop,2,mean)
lines(mean.traj, lwd=2)

# save the outputs to file. Note this currently gets overwritten each time the script runs.

save(all.outputs, file='jpopmodel_data.Rdata')

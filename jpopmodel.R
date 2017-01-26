# R script to implement a PVA model for Brazilian Jaguar in the Atlandtic Forest region.
# Author: Ascelin Gordon ascelin.gordon@rmit.edu.au
# Commencement Date: 02 Feb 2016


# To run:
#   source( 'jpopmodel.R')

rm(list=ls(all=TRUE))

set.seed(1)

source('jpopmodel.functions.R')


        # ------------------------------------------------
        # Define the name parameters of the model
        # ------------------------------------------------

DEBUG.LEVEL <- 1  # 0-none; 1-terse, 2-verbose
num.time.steps <- 10
num.reps <- 50
include.dispersal <- FALSE
output.filename <- 'jpopmodel_data_nodisp.Rdata'


        # ------------------------------------------------
        # It's planned that These parameters will be eventually be
        # extracted from the GIS and other data data we will be
        # given. But just specifying them here for now
        # ------------------------------------------------


num.jcus <- 4
num.life.stages <- 3

jcu.att <- data.frame( cc=c(5,10,7,15),
                       #cc=c(5,45,45,40),
                      #mortality=c(0.5,0.15,0.2,0.4), #assuming same for each stage for now so one no per JCU
                      mortality=rep(0.2,num.jcus),
                      
                      # Assuming only last stage gives birth, this is the prob of giving birth to 1 offspring
                      birth.rate=rep(0.2, num.jcus)
                      #birth.rate=c(0.1, 0.2, 0.25, 0.28)
                      )

# Make an initial population (random from now) 
initial.pop <- matrix( ncol=num.jcus, nrow=num.life.stages, sample(1:10, size=num.jcus*3, replace=TRUE ) )
#initial.pop <- matrix( ncol=num.jcus, nrow=num.life.stages, c(10,7,10,12),byrow=TRUE )
colnames(initial.pop) <- paste( 'jcu', 1:num.jcus, sep='') # set column names
rownames(initial.pop) <- paste( 'stage', 1:num.life.stages, sep='') # set the row names

# Make a dispersal matrxi between each of the JCUs
disp.mort.mat <- matrix(ncol=num.jcus, nrow=num.jcus)
disp.mort.mat[] <- round(runif(min=0.2, max=0.6, num.jcus^2),2)  # set all values to be the same for now
diag(disp.mort.mat) <- 0  # set the diagonal values to zero.


        # ------------------------------------------------
        # Validity checks
        # ------------------------------------------------

if( dim(jcu.att)[1] != num.jcus ) stop( '\n\nERROR JCU attributes not conistent with number of JCUs')


        # ------------------------------------------------
        # Make some objects to store results as the model runs
        # ------------------------------------------------

# A data.frame to store the full state of the system
all.outputs <- data.frame ( "rep"=NA, "time"=NA, "jcu"=NA, "cc"=NA, "mortality"=NA, "birth.rate"=NA,
                           "stage1"=NA, "stage2"=NA, "stage3"=NA ) [numeric(0), ]

# A datafram to just store the total population size to make a quick plot at the end
total.pop <- matrix( ncol = num.time.steps, nrow = num.reps )

    
    
        # ------------------------------------------------
        # Run the model
        # ------------------------------------------------


for( rep in 1:num.reps) {
    
    if(DEBUG.LEVEL>0){cat( "\n\n--------------------------------")}
    cat( '\n Rep = ', rep )
    
    # First save the initial values for timestep 1 for all jcus
    time<-1
    for( jcu in 1:num.jcus){
        init.info <- c(rep, time, jcu, jcu.att[jcu,"cc"], jcu.att[jcu,"mortality"], 
                  jcu.att[jcu,"birth.rate"], initial.pop[,jcu] )

        # If first data entry, replace the first line
        if( dim(all.outputs)[1] == 0 ) all.outputs[1,] <- init.info
        else  all.outputs <- rbind( all.outputs, init.info ) # otherwise append it           
    }
    # Save the total population size to make a quick plot at the end 
    total.pop[rep, time] <- sum(initial.pop)
    if(DEBUG.LEVEL>0) cat( '\n Time step = ', time )
    if(DEBUG.LEVEL>1) show(initial.pop) 
   

    # Now loop over time steps from 2 onwards
    current.pop <- initial.pop
    for( time in 2:num.time.steps) {
        

        if(DEBUG.LEVEL>0) cat( '\n Time step = ', time)

        # age the population
        current.pop <- age.population(current.pop, num.life.stages)
            
        # do reproduction (assumes for now only last stage reproduces)
        current.pop <- reproduce(current.pop, num.life.stages, jcu.att$birth.rate)
        
        # loop through each JCU and apply the JCU specific mortality to each life stage
        for( jcu in 1:num.jcus){

            # Apply mortality to the populations in each JCU
            current.pop[,jcu] <- apply.mortality(current.pop[,jcu], jcu.att$mortality[jcu])
        }

        # apply dispersal
        if( include.dispersal ) 
            current.pop <- apply.dispersal(current.pop, jcu.att$cc, disp.mort.mat)
        
        # ------------------------------------------------
        # Save the current info 
        # ------------------------------------------------

        # Save the total population size to make a quick plot at the end
        total.pop[rep, time] <- sum(current.pop)

        for( jcu in 1:num.jcus){        
            # Build a vector for the current info of the system
            current.info <- c(rep, time, jcu, jcu.att[jcu,"cc"], jcu.att[jcu,"mortality"], 
                	           jcu.att[jcu,"birth.rate"], current.pop[,jcu] )
               	
            # Add the the all outouts dataframe
            all.outputs <- rbind( all.outputs, current.info )
        }

        #if(DEBUG){ show(current.pop) }
        
        
        
    }
    
    if(DEBUG.LEVEL>0) cat("\n")
  
}


        # ------------------------------------------------
        # Make some simple plots. For more detailed analysis run the
        # analyse.results.R script
        # ------------------------------------------------


# Plot the total pop trajectory of each realisation and also the mean
# trajectory (using matplot() to automatically plot a curve for each
# realization)
matplot ( t(total.pop), type = 'l', main= 'Total population', xlab='time', ylab='pop size',
         ylim=c(0, max(total.pop)) )
mean.traj <- apply(total.pop,2,mean)
lines(mean.traj, lwd=2)

# Save the outputs to file. This is used in the analyse.results.R
# script. Note this currently gets overwritten each time the script
# runs.
save(all.outputs, file=output.filename)

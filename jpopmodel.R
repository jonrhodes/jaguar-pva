# R script to implement a PVA model for Brazilian Jaguar in the Atlandtic Forest region.
# Author: Ascelin Gordon ascelin.gordon@rmit.edu.au
# Commencement Date: 02 Feb 2016


# To run:
#   source( 'jpopmodel.R')

rm(list=ls(all=TRUE))

set.seed(2)


source('jpopmodel.functions.R')



        # ------------------------------------------------
        # Define the name parameters of the model
        # ------------------------------------------------

DEBUG.LEVEL <- 1  # 0-none; 1-terse, 2-verbose
# num.time.steps <- 80
# num.reps <- 20
num.time.steps <- 20
num.reps <- 1
OPT.INCLUDE.DISPERSAL <- TRUE
output.filename <- 'jpopmodel_data_disp.Rdata'

# Option of whether to use a truncated Poisson distribution for determing the
# number offspring (TRUE) or whether to assume there are only 1 or 2 offspring
# (FALSE). The truncated Poisson means that occasionally there may be 5 or 6
# offspring.
OPT.USE.TRUNCATED.POISSON.FOR.REPRODUCTION <- TRUE 


        # ------------------------------------------------
        # It's planned that These parameters will be eventually be
        # extracted from the GIS and other data data we will be
        # given. But just specifying them here for now
        # ------------------------------------------------


num.jcus <- 3
num.life.stages <- 3


jcu.att <- data.frame( #cc=c(20,20,20), #15),
                       cc=c(4,4,4),

                      #mortality.stage3=c(0.5,0.15,0.2,0.4), #assuming same for each stage for now so one no per JCU
                      # mortality.stage1=rep(0.38,num.jcus),
                      # mortality.stage2=rep(0.26,num.jcus),
                      # mortality.stage3=rep(0.14,num.jcus),
                      mortality.stage1=c(0.38, 0.44, 0.58),
                      mortality.stage2=c(0.26, 0.33, 0.55),
                      mortality.stage3=c(0.14, 0.27, 0.51),
                      mortality.floaters=c(0.14, 0.27, 0.51), # assume floaters the same as stage 3 for now.

                      # Assuming only last stage gives birth, this is the prob of giving birth to 1 or more offspring
                      birth.rate.mean=rep(0.45, num.jcus), 
                      birth.rate.upper.bound=rep(0.3, num.jcus), # not currently used
                      birth.rate.lower.bound=rep(0.6, num.jcus)  # not currently used
                      
                      )

# Make an initial population (random from now), row for each 3 life stage and one for the floaters 
#initial.pop <- matrix( ncol=num.jcus, nrow=num.life.stages, sample(1:10, size=num.jcus*3, replace=TRUE ) )
initial.pop <- matrix( ncol=num.jcus, nrow=num.life.stages+1, rep(5, num.jcus*4) ) # in in each stage

#initial.pop <- matrix( ncol=num.jcus, nrow=num.life.stages, c(10,7,10,12),byrow=TRUE )
colnames(initial.pop) <- paste0( 'jcu', 1:num.jcus) # set column names
rownames(initial.pop) <- c(paste0( 'stage', 1:num.life.stages ), 'floaters' ) # set the row names

# Set the number of floaters to be zero
initial.pop['floaters',] <- 0



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
all.outputs <- data.frame ( "rep"=NA, "time"=NA, "jcu"=NA, "cc"=NA, "mortality.stage3"=NA, "birth.rate.mean"=NA,
                           "stage1"=NA, "stage2"=NA, "stage3"=NA, "floaters"=NA ) [numeric(0), ]

# A dataframx to just store the total population size to make a quick plot at the end
total.pop <- matrix( ncol = num.time.steps, nrow = num.reps )

    
    
        # ------------------------------------------------
        # Run the model
        # ------------------------------------------------


for( rep in 1:num.reps) {
    
    if(DEBUG.LEVEL>0){cat( "\n\n--------------------------------")}
    if(rep%%10==0) cat( '\n Rep = ', rep ) # print every 10th rep 
    
    # First save the initial values for timestep 1 for all jcus
    time<-1
    for( jcu in 1:num.jcus){
        init.info <- c(rep, time, jcu, jcu.att[jcu,"cc"], jcu.att[jcu,"mortality.stage3"], 
                  jcu.att[jcu,"birth.rate.mean"], initial.pop[,jcu] )

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
        current.pop <- reproduce(current.pop, jcu.att$birth.rate.mean, litter.size.dist)
        
        # loop through each JCU and apply the JCU specific mortality to each life stage
        for( jcu in 1:num.jcus){
            #browser()
            # Apply mortality to the populations in each JCU
            current.pop[,jcu] <- apply.mortality(current.pop[,jcu], jcu.att$mortality.stage1[jcu],
                                                jcu.att$mortality.stage2[jcu],
                                                jcu.att$mortality.stage3[jcu], 
                                                jcu.att$mortality.floaters[jcu])
        }

        # apply dispersal
        if( OPT.INCLUDE.DISPERSAL ) 
            current.pop <- apply.dispersal(current.pop, jcu.att$cc, disp.mort.mat)
        
        # ------------------------------------------------
        # Save the current info 
        # ------------------------------------------------

        # Save the total population size to make a quick plot at the end
        total.pop[rep, time] <- sum(current.pop)

        for( jcu in 1:num.jcus){        
            # Build a vector for the current info of the system
            current.info <- c(rep, time, jcu, jcu.att[jcu,"cc"], jcu.att[jcu,"mortality.stage3"], 
                	           jcu.att[jcu,"birth.rate.mean"], current.pop[,jcu] )
               	
            # Add the the all outouts dataframe
            all.outputs <- rbind( all.outputs, current.info )
        }

        #if(DEBUG){ show(current.pop) }
        
        
      if(DEBUG.LEVEL>0) cat("\n-----------------------------------------------\n")    
    }
    
    
  
}
 
cat('\n')


        # ------------------------------------------------
        # Make some simple plots. For more detailed analysis run the
        # analyse.results.R script
        # ------------------------------------------------


# Plot the total pop trajectory of each realisation and also the mean
# trajectory (using matplot() to automatically plot a curve for each
# realization)
matplot ( t(total.pop), type = 'l', main= 'Total population', xlab='time', ylab='pop size',
         #ylim=c(0, 120)
         ylim=c(0, max(total.pop))
        )
mean.traj <- apply(total.pop,2,mean)
lines(mean.traj, lwd=2)

# Save the outputs to file. This is used in the analyse.results.R
# script. Note this currently gets overwritten each time the script
# runs.
save(all.outputs, file=output.filename)

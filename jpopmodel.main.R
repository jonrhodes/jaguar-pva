

run.jpop.model <-function(initial.pop, jcu.att, disp.mort.mat) {


            # ------------------------------------------------
            # Validity checks
            # ------------------------------------------------

    if( dim(jcu.att)[1] != num.jcus ) stop( '\n\nERROR JCU attributes not conistent with number of JCUs')


            # ------------------------------------------------
            # Make some objects to store results as the model runs
            # ------------------------------------------------

    # A data.frame to store the full state of the system
    all.outputs <- data.frame ( "rep"=NA, "time"=NA, "jcu"=NA, "K"=NA, "mortality.stage3"=NA, "birth.rate.mean"=NA,
                               "stage1"=NA, "stage2"=NA, "stage3"=NA, "floaters"=NA ) [numeric(0), ]

    # A dataframx to just store the total population size to make a quick plot at the end
    total.pop <- matrix( ncol = num.time.steps, nrow = num.reps )

        
        
            # ------------------------------------------------
            # Run the model
            # ------------------------------------------------


    for( rep in 1:num.reps) {
        
        if(DEBUG.LEVEL>1){cat( "\n\n--------------------------------")}
        if(rep%%10==0) cat( '\n Rep = ', rep ) # print every 10th rep 
        
        # First save the initial values for timestep 1 for all jcus
        time<-1
        for( jcu in 1:num.jcus){
            init.info <- c(rep, time, jcu, jcu.att[jcu,"K"], jcu.att[jcu,"mortality.stage3"], 
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
                current.pop <- apply.dispersal(current.pop, jcu.att$K, disp.mort.mat)
            
            # ------------------------------------------------
            # Save the current info 
            # ------------------------------------------------

            # Save the total population size to make a quick plot at the end
            total.pop[rep, time] <- sum(current.pop)

            for( jcu in 1:num.jcus){        
                # Build a vector for the current info of the system
                current.info <- c(rep, time, jcu, jcu.att[jcu,"K"], jcu.att[jcu,"mortality.stage3"], 
                                   jcu.att[jcu,"birth.rate.mean"], current.pop[,jcu] )
                    
                # Add the the all outouts dataframe
                all.outputs <- rbind( all.outputs, current.info )
            }

            #if(DEBUG){ show(current.pop) }
            
            
          if(DEBUG.LEVEL>1) cat("\n-----------------------------------------------\n")    
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
}


            # ------------------------------------------------
            # Low level options
            # ------------------------------------------------

DEBUG.LEVEL <- 1  # 0-none; 1-terse, 2-verbose

# Option of whether to use a truncated Poisson distribution for determining the
# number offspring (TRUE) or whether to assume there are only 1 or 2 offspring
# (FALSE). The truncated Poisson means that occasionally there may be 5 or 6
# offspring.
OPT.USE.TRUNCATED.POISSON.FOR.REPRODUCTION <- TRUE
OPT.INCLUDE.DISPERSAL <- TRUE
OPT.NUMBER.OF.LIFE.STAGES <- 3

run.jpop.model <-function(expert.ID, expert.realization, num.stoch.realizatons, initial.pop, jcu.att, disp.mort.mat, num.time.steps) {

    # determine the number of JCUs from the jcu.attributes data frame
    num.jcus <- dim(jcu.att)[1]

            # ------------------------------------------------
            # Make some objects to store results as the model runs
            # ------------------------------------------------

    # A data.frame to store the full state of the system
    all.outputs <- data.frame ( 'expert.ID'=NA, 'expert.realization'=NA, 'stoch.realization'=NA, 'time'=NA, 'jcu'=NA, 'K'=NA, 'mortality.stage3'=NA, 'birth.rate.mean'=NA,
                               'stage1'=NA, 'stage2'=NA, 'stage3'=NA, 'floaters'=NA ) [numeric(0), ]

    # A dataframx to just store the total population size to make a quick plot at the end
    total.pop <- matrix( ncol = num.time.steps, nrow = num.stoch.realizatons )



            # ------------------------------------------------
            # Run the model
            # ------------------------------------------------


    for( rep in 1:num.stoch.realizatons) {

        if(DEBUG.LEVEL>1){cat( "\n\n--------------------------------")}
        if(rep%%10==0) cat( '\n Rep = ', rep ) # print every 10th rep

        # First save the initial values for timestep 1 for all jcus
        time<-1
        for( jcu in 1:num.jcus){
            init.info <- c(expert.ID, expert.realization, rep, time, jcu, jcu.att[jcu,"K"], jcu.att[jcu,"mortality.stage3"],
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
            current.pop <- age.population(current.pop, OPT.NUMBER.OF.LIFE.STAGES)

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
                current.info <- c(expert.ID, expert.realization, rep, time, jcu, jcu.att[jcu,"K"], jcu.att[jcu,"mortality.stage3"],
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
    plot.title <- paste0('Total population (expert ', expert.ID, ')' )
    matplot ( t(total.pop), type = 'l', main= plot.title, xlab='time', ylab='pop size',
             #ylim=c(0, 120)
             ylim=c(0, max(total.pop))
            )
    mean.traj <- apply(total.pop,2,mean)
    lines(mean.traj, lwd=2)

    # Save the outputs to file. This is used in the analyse.results.R
    # script. Note this currently gets overwritten each time the script
    # runs.

    #save(all.outputs, file=output.filename)

    return (all.outputs)
}

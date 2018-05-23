

# To run
#    source( 'analyse.results.R' )

rm(list=ls(all=TRUE))


# read in the results to analyse 
results.file <- 'jpopmodel_data_disp.Rdata'
#results.file <- 'jpopmodel_test.Rdata'
model.output.all.experts <- readRDS ( results.file )


# as have added data for running the model with multiple experts, just plot
# the results for expert 1 for now. Later will want to generate plots for each
# expert
# cat('\nPloting results only for exper 1 for now')
# model.output <- subset( model.output, expert.ID==1)


expert.vec <- unique(model.output.all.experts$expert.ID)
expert.realization.vec <- unique(model.output.all.experts$expert.realization)

cat('Ploting resuts for expert for', length(expert.vec), 'experts\n' )


for(current.expert in expert.vec) {

    cat('\nPloting results only for expert ID', current.expert)

    # pull out the result for the given expert
    model.output <- subset( model.output.all.experts, expert.ID==current.expert)

    

    # Get vales from the dataframe
    num.reps <- max(model.output$stoch.realization)
    time.steps.vec <- unique(model.output$time)
    max.time <- max(time.steps.vec)
    jcu.vec <- unique(model.output$jcu)



            # ------------------------------------------------
            # Analyse the resuls
            # ------------------------------------------------


    # Extract out and plot the total population change over time

    # Make a matrix to hold the data
    pop.traj <- matrix(ncol=length(time.steps.vec), nrow=num.reps)


    # Made an array to hold the data for all JCUs
    pop.traj.jcu.array <- array(dim=c(length(time.steps.vec), nrow=num.reps, length(jcu.vec) ) )

    # Below is extracting the data for each JCU, could do the same for each
    # life stage also

    cat( '\nExtracting data for ecach JCU...')

    # Extract out the data to plot
    for(i in 1:num.reps ) {
        
        if(i%%10==0) cat( '\n Rep = ', i) # print every 10th rep 

        for( x in time.steps.vec ){

            pop.traj[i, x] <- sum( subset( model.output, stoch.realization==i & time==x, select=total.pop ) )

            # Pull out the results for each JCU
            for( cur.jcu in jcu.vec ) {

                pop.traj.jcu.array[x,i,cur.jcu] <- sum( subset( model.output, stoch.realization==i & time==x & jcu==cur.jcu, select=total.pop ) )

            }

        }
    }

    cat('\nGenerating plots')

    par(mfrow=c(3,3))

    tot.cc <- sum(subset( model.output, stoch.realization==1 & time==1, select=K ))

    # Use matplot to plot a curve for each realization 
    mean.traj <- apply(pop.traj,2,mean)
    main.txt <- paste('Tot pop (MEP=', mean.traj[max.time],')', sep='')
    matplot(t(pop.traj), type='l', main=main.txt, xlab='time', ylab='pop size',
        ylim=c(0, max(tot.cc, pop.traj)) )
    plot.title <- paste0(results.file, ' (expert ', current.expert, ')')
    mtext(plot.title, outer=TRUE, line=-1.5)

    # plot the mean trajectory
    lines(mean.traj, lwd=3)

    # add a line for the total carrying capacity
    #abline(h=tot.cc, col='grey')

    # The abline call justs plots the carrying capacity for the jcu (just the
    # initial value for now)

    # jcus.to.plot <- jcu.vec
    jcus.to.plot <- 5:17
    jcus.to.plot <- c(5,9,14,15)

    for(cur.jcu in jcus.to.plot ) {

        # plot jcu 1
        mean.traj <- apply(pop.traj.jcu.array[,,cur.jcu],1,mean)

        #main.txt <- paste('JCU1 (MEP=', mean.traj[max.time],')', sep='')
        main.txt <- paste0('JCU ', cur.jcu) 
        matplot(pop.traj.jcu.array[,,cur.jcu], type='l', main=main.txt, xlab='time', ylab='pop size' )
        lines(mean.traj, lwd=3)
        #abline( h=subset( model.output, jcu==1 & stoch.realization==1 & time==1, select=cc ), col='grey' )
        

    }

    cat('\n')
    

    

}

par(mfrow=c(1,1))



# To run
#    source( 'analyse.results.R' )

rm(list=ls(all=TRUE))


# Creates dataframe "model.output"
results.file <- 'jpopmodel_data_disp.Rdata'
load ( results.file )


# Get vales from the dataframe
num.reps <- max(model.output$stoch.realization)
time.steps.vec <- unique(model.output$time)
max.time <- max(time.steps.vec)


        # ------------------------------------------------
        # Analyse the resuls
        # ------------------------------------------------


# Extract out and plot the total population change over time

# Make a matrix to hold the data
pop.traj <- matrix(ncol=length(time.steps.vec), nrow=num.reps)

# Made a matrix to hold the data for the first 3 JCUs
pop.traj.jcu1 <- matrix(ncol=length(time.steps.vec), nrow=num.reps) 
pop.traj.jcu2 <- matrix(ncol=length(time.steps.vec), nrow=num.reps) 
pop.traj.jcu3 <- matrix(ncol=length(time.steps.vec), nrow=num.reps) 

# Below is extracting the data for each JCU, could do the same for each
# life stage also

# Extract out the data to plot
for(i in 1:num.reps ) {
    for( x in time.steps.vec ){
        tmp <- subset( model.output, stoch.realization==i & time==x, select=stage1:stage3 )      
        pop.traj[i, x ] <- sum(tmp)

        tmp1 <- subset( model.output, stoch.realization==i & time==x & jcu==1, select=stage1:stage3 )      
        pop.traj.jcu1[i, x ] <- sum(tmp1)
        
        tmp2 <- subset( model.output, stoch.realization==i & time==x & jcu==2, select=stage1:stage3 )      
        pop.traj.jcu2[i, x ] <- sum(tmp2)
        
        tmp3 <- subset( model.output, stoch.realization==i & time==x & jcu==3, select=stage1:stage3 )      
        pop.traj.jcu3[i, x ] <- sum(tmp3)

    }
}

par(mfrow=c(2,2))

tot.cc <- sum(subset( model.output, stoch.realization==1 & time==1, select=K ))

# Use matplot to plot a curve for each realization 
mean.traj <- apply(pop.traj,2,mean)
main.txt <- paste('Tot pop (MEP=', mean.traj[max.time],')', sep='')
matplot(t(pop.traj), type='l', main=main.txt, xlab='time', ylab='pop size',
    ylim=c(0, max(tot.cc, pop.traj)) )
mtext(results.file, outer=TRUE, line=-1.5)

# plot the mean trajectory
lines(mean.traj, lwd=3)

# add a line for the total carrying capacity
#abline(h=tot.cc, col='grey')

# The abline call justs plots the carrying capacity for the jcu (just the
# initial value for now)

# plot jcu 1
mean.traj <- apply(pop.traj.jcu1,2,mean)
#main.txt <- paste('JCU1 (MEP=', mean.traj[max.time],')', sep='')
main.txt <- paste('JCU1 Low hunting pressure') 
matplot(t(pop.traj.jcu1), type='l', main=main.txt, xlab='time', ylab='pop size' )
lines(mean.traj, lwd=3)
#abline( h=subset( model.output, jcu==1 & stoch.realization==1 & time==1, select=cc ), col='grey' )

# plot jcu 2
mean.traj <- apply(pop.traj.jcu2,2,mean)
#main.txt <- paste('JCU2 (MEP=', mean.traj[max.time],')', sep='')
main.txt <- paste('JCU2 Medium hunting pressure')
matplot(t(pop.traj.jcu2), type='l', main=main.txt, xlab='time', ylab='pop size')
lines(mean.traj, lwd=3)
#abline( h=subset( model.output, jcu==2 & stoch.realization==1 & time==1, select=cc ), col='grey' )

# plot jcu 3
mean.traj <- apply(pop.traj.jcu3,2,mean)
#main.txt <- paste('JCU3 (MEP=', mean.traj[max.time],')', sep='')
main.txt <- paste('JCU3 High hunting pressure')
matplot(t(pop.traj.jcu3), type='l', main=main.txt, xlab='time', ylab='pop size')
lines(mean.traj, lwd=3)
#abline( h=subset( model.output, jcu==3 & stoch.realization==1 & time==1, select=cc ), col='grey' )


par(mfrow=c(1,1))


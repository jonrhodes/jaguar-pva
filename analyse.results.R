

# To run
#    source( 'analyse.results.R' )

rm(list=ls(all=TRUE))


# creates dataframe "all.outputs"
load ( 'jpopmodel_data.Rdata' )


# get vales from the dataframe
num.reps <- max(all.outputs$rep)
time.steps.vec <- unique(all.outputs$time)



        # ------------------------------------------------
        # Analyse the resuls
        # 
        # ------------------------------------------------


# extract out and plot the total population change over time

# make a matrix to hold the data
pop.traj <- matrix(ncol=length(time.steps.vec), nrow=num.reps)

# made a matrix to hold the data for the first 3 JCUs
pop.traj.jcu1 <- matrix(ncol=length(time.steps.vec), nrow=num.reps) 
pop.traj.jcu2 <- matrix(ncol=length(time.steps.vec), nrow=num.reps) 
pop.traj.jcu3 <- matrix(ncol=length(time.steps.vec), nrow=num.reps) 

# extract out the 
for(i in 1:num.reps ) {
    for( x in time.steps.vec ){
        tmp <- subset( all.outputs, rep==i & time==x, select=stage1:stage3 )      
        pop.traj[i, x ] <- sum(tmp)

        tmp1 <- subset( all.outputs, rep==i & time==x & jcu==1, select=stage1:stage3 )      
        pop.traj.jcu1[i, x ] <- sum(tmp1)
        
        tmp2 <- subset( all.outputs, rep==i & time==x & jcu==2, select=stage1:stage3 )      
        pop.traj.jcu2[i, x ] <- sum(tmp2)
        
        tmp3 <- subset( all.outputs, rep==i & time==x & jcu==3, select=stage1:stage3 )      
        pop.traj.jcu3[i, x ] <- sum(tmp3)
        
    }
}

par(mfrow=c(2,2))

# plot all curves
matplot(t(pop.traj), type="l")

# plot the mean trajectory
mean.traj <- apply(pop.traj,2,mean)
lines(mean.traj, lwd=3)


# plot jcu 1
matplot(t(pop.traj.jcu1), type="l")
mean.traj <- apply(pop.traj.jcu1,2,mean)
lines(mean.traj, lwd=3)

# plot jcu 2
matplot(t(pop.traj.jcu2), type="l")
mean.traj <- apply(pop.traj.jcu2,2,mean)
lines(mean.traj, lwd=3)

# plot jcu 3
matplot(t(pop.traj.jcu3), type="l")
mean.traj <- apply(pop.traj.jcu3,2,mean)
lines(mean.traj, lwd=3)

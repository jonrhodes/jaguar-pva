sim.num.JCUs <- 3
num.life.stages <- 3

jcu.attributes <- data.frame( K=c(20,20,20), #15),
                       #K=c(4,4,4),

                      #mortality.stage3=c(0.5,0.15,0.2,0.4), #assuming same for each stage for now so one no per JCU
                      # mortality.stage1=rep(0.38,sim.num.JCUs),
                      # mortality.stage2=rep(0.26,sim.num.JCUs),
                      # mortality.stage3=rep(0.14,sim.num.JCUs),
                      mortality.stage1=c(0.38, 0.44, 0.58),
                      mortality.stage2=c(0.26, 0.33, 0.55),
                      mortality.stage3=c(0.14, 0.27, 0.51),
                      mortality.floaters=c(0.14, 0.27, 0.51), # assume floaters the same as stage 3 for now.


                      birth.rate.mean=rep(0.45, sim.num.JCUs)
                      )


# Make an initial population (random from now), row for each 3 life stage and one for the floaters
#initial.population <- matrix( ncol=sim.num.JCUs, nrow=num.life.stages, sample(1:10, size=sim.num.JCUs*3, replace=TRUE ) )
initial.population <- matrix( ncol=sim.num.JCUs, nrow=num.life.stages+1, rep(20, sim.num.JCUs*4) ) # in in each stage

#initial.population <- matrix( ncol=sim.num.JCUs, nrow=num.life.stages, c(10,7,10,12),byrow=TRUE )
colnames(initial.population) <- paste0( 'jcu', 1:sim.num.JCUs) # set column names
rownames(initial.population) <- c(paste0( 'stage', 1:num.life.stages ), 'floaters' ) # set the row names

# Set the number of floaters to be zero
initial.population['floaters',] <- 0

# Make a dispersal matrxi between each of the JCUs
disp.mort.matrix <- matrix(ncol=sim.num.JCUs, nrow=sim.num.JCUs)
disp.mort.matrix[] <- round(runif(min=0.2, max=0.6, sim.num.JCUs^2),3)  # set all values to be the same for now
diag(disp.mort.matrix) <- 0  # set the diagonal values to zero.

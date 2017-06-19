
#---------------------------------------------------------------

rtpois <- function( num.of.samples, pre.truncated.mean ) {

	# Function to sample from a truncated Poisson distribution (ie with
	# excluding zero value). I got his code from the following refrence:
	# http://tolstoy.newcastle.edu.au/R/help/05/05/3745.html  Note that given
	# a pre-trucated mean of T, then the mean of the truncated distribution
	# will be given by T/(1 - exp(-T)). For the litter size, the expert
	# elicited values are for the mean of the trucated distibution (eg mean
	# litter size of 1.62 means for those give birth, the mean number of
	# offspring will be 1.62). As this equation is hard to solve for T, I've
	# made a look up table here:
	#   pre.truncated.mean of 1.06 gives a truncated mean of 1.62

	T <- pre.truncated.mean
  	U<-runif(num.of.samples)   # the uniform sample

	if(DEBUG.LEVEL>1){	cat('\n The resulting trucated mean is', T/(1 - exp(-T)), '\n' ) }

  	t = -log(1 - U*(1 - exp(-T))) # the "first" event-times   
  	T1<-(T - t) # the set of (T-t)
  	X <- rpois(num.of.samples,T1)+1  # the final truncated Poisson sample
  	return(X)

}

#---------------------------------------------------------------

age.population <- function( cur.pop, num.stages) {

    # Move the population through the age stages

    # Note: curr.pop is matrix of life stages (rows) by JCUs (columns), showing
    # the number of individuals in each life stage in each JCU.
    
    age = function(x) {
        tmp <- x
        for( i in num.stages:1) {
            if( i == num.stages) x[i] <- x[i] + x[i-1] # last stage
            else if (i == 1)  x[i] <- 0 # first stage
            else x[i] <- x[i-1] # other stages  
        }         
        return (x)
    }
    
    return ( apply(cur.pop,2,age)  )
}

#---------------------------------------------------------------

reproduce <- function( cur.pop, birth.rate, litter.size.dist ) {

    
    s1 <- cur.pop[1,]
    s3 <- cur.pop[3,]
    pop <- rep(0, length(s1) )

    # Loop over JCUs and decide which give birth
    for( i in 1:length(s1) ) {

    
    	# First determine the number in stage 3 that reproduce
    	num.that.reproduce <- rbinom(n=1, size=s3[i], prob=birth.rate[i])

    	# Now determine how many offspring each one has
    	
    	if( OPT.USE.TRUNCATED.POISSON.FOR.REPRODUCTION) {
    		# Use the truncated Poisson distribution (defined above). Note that
    		# pre.truncated.mean or 1.06 correcponds to runcated.mean of 1.62
    		num.new.offspring <- sum(rtpois(num.that.reproduce, pre.truncated.mean=1.06))

    	} else {
    		# In this case assume that jaguar litter size is only 1 or 2
    		# offspring. Set the probabilities of having 1 or 2 offspring such
    		# that expected litter size has a value of 1.62 (the mean expert
    		# estimate), which then generates 0.38 * 1 + (1-0.38)*2 = 1.62
			litter.size.dist <- c(prob.litter.size.eq.1=0.38, prob.litter.size.eq.2=1-0.38)

    		num.new.offspring <- sum(sample(c(1,2), size=num.that.reproduce, replace=TRUE, prob=litter.size.dist))

    	}
    	#browser()
    
    	pop[i] <- s1[i] + num.new.offspring	

    }

#browser()
    cur.pop[1,] <- pop
    
    return ( cur.pop  )
}

#---------------------------------------------------------------
apply.mortality <- function(pop, mortality.s1, mortality.s2, mortality.s3 ){

    # flip a weighted coin to decide how many die
    # note: currently assumes all stages have the same mortality
    # mortality.function = function(x) rbinom(n=1, size=x, prob=(1-mortality))
    # return( sapply( pop, mortality.function) )

    pop['stage1'] <- rbinom(n=1, size=pop['stage1'], prob=(1-mortality.s1))
    pop['stage2'] <- rbinom(n=1, size=pop['stage2'], prob=(1-mortality.s2))    
    pop['stage3'] <- rbinom(n=1, size=pop['stage3'], prob=(1-mortality.s3))

    return(pop)
    
}

#---------------------------------------------------------------

disperse.stage <- function() {


}

#---------------------------------------------------------------

apply.dispersal <- function( current.pop, jcu.cc, disp.mort.mat ) {

	# For now, all jaguars over cc try and disperse

	if(DEBUG.LEVEL > 1) {cat( '\n  In dispersal function\n'); show(current.pop) }

	jcu.vec <- 1:length(jcu.cc)
	dispersal.ctr <- 0
	dispersal.mort.ctr <- 0

	for( i in jcu.vec ){

		if(DEBUG.LEVEL > 1) cat('\njcu=', i, 'cc=', jcu.cc[i], 'stage3.pop=',current.pop[3,i], 
		 	'no over cc=', max(current.pop[3,i] - jcu.cc[i],0) , '\n')

		source.jcu <- i
		num.to.disperse <- current.pop[3,i] - jcu.cc[i]

		
		# If there are more stage 3 adults than the cc disperse them
		if( num.to.disperse > 0) {


			# choose the jcus that each jaguar try and disperse to 
			# select each one randomly for now
			dest.jcus <- sample(jcu.vec[-i], num.to.disperse, replace=TRUE )

			# reduce the pop of the source jcu by the number that disperse
			current.pop[3,i] <- current.pop[3,i] - num.to.disperse

			# determine which ones survive the dispersal
			survival.vec <- rbinom(n=num.to.disperse, size=1, prob=(1-disp.mort.mat[source.jcu, dest.jcus]) )
			dest.jcus.surviving <- dest.jcus[which(survival.vec==1)]
			no.die.dispersing <- length(which(survival.vec==0))

			# determine the total number going to each jcu
			# Note: rle: Run Length Encoding to compute the lengths and values of runs of equal values 
			#            in a vector
			rle.cts <- rle(sort(dest.jcus.surviving))
			receive.cts <- rep(0, length(jcu.vec))
			receive.cts[rle.cts$values] <-  rle.cts$lengths

			# Update the population
			current.pop[3,] <- current.pop[3,] + receive.cts

			# track some dispersal stats
			dispersal.ctr <- dispersal.ctr + num.to.disperse
			dispersal.mort.ctr <- dispersal.mort.ctr + no.die.dispersing

			
		}


	}
	#browser()
	
	if(DEBUG.LEVEL>0) cat(' (No disp:', dispersal.ctr, 'mort in disp:', dispersal.mort.ctr, ')')
	
	return (current.pop)
}

#---------------------------------------------------------------

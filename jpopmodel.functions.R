
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
  	U <- runif(num.of.samples)   # the uniform sample

	#if(DEBUG.LEVEL>1){	cat('\n The resulting trucated mean is', T/(1 - exp(-T))) }

  	t = -log(1 - U*(1 - exp(-T))) # the "first" event-times   
  	T1 <- (T - t) # the set of (T-t)
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

    pop <- rep(0, length(s3) )

    # Loop over JCUs and decide which give birth
    for( i in 1:length(s3) ) {

    
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

apply.mortality <- function(pop, mortality.s1, mortality.s2, mortality.s3, 
                            mortality.floaters ){

    # flip a weighted coin to decide how many die
    #browser()
    pop['stage1'] <- rbinom(n=1, size=pop['stage1'], prob=(1-mortality.s1))
    pop['stage2'] <- rbinom(n=1, size=pop['stage2'], prob=(1-mortality.s2))    
    pop['stage3'] <- rbinom(n=1, size=pop['stage3'], prob=(1-mortality.s3))
    pop['floaters'] <- rbinom(n=1, size=pop['floaters'], prob=(1-mortality.floaters))

    #browser()
    return(pop)
    
}

#---------------------------------------------------------------
# function to work out how many floaters can disperse into their 
# attached JCU
determine.floaters.moving.to.jcu <- function(K, pop, floaters){

    territories.available <- max(K - pop, 0)
    
    num.floaters.moving.into.jcu <- 0

    if( floaters > 0 & territories.available > 0){
            # if that is the case, floaters first disperse into available
            # territories of the JCU they are attached to
            if( territories.available >= floaters ){    
                # all floaters move inot the JCU
                num.floaters.moving.into.jcu <- floaters
            } else {
                # there are more floaters than territories
                num.floaters.moving.into.jcu <- territories.available
            }
    }


        return(num.floaters.moving.into.jcu)
}

#tests

#determine.floaters.moving.to.jcu(K=4, pop=2, floaters=4 ) # ans=2
#determine.floaters.moving.to.jcu(K=4, pop=4, floaters=4 ) # ans=0
#determine.floaters.moving.to.jcu(K=10, pop=15, floaters=15 ) # ans 0
#determine.floaters.moving.to.jcu(K=10, pop=0, floaters=15 ) # ans 10

#---------------------------------------------------------------

apply.dispersal <- function( pop, jcu.cc, disp.mort.mat ) {

	if(DEBUG.LEVEL > 1) {cat( '\n  In dispersal function\n'); show(pop) }

	jcu.vec <- 1:length(jcu.cc)
	dispersal.ctr <- 0
	dispersal.mort.ctr <- 0

	for( source.jcu in jcu.vec ){

        # print out current JCU  
        i<- source.jcu
        if(DEBUG.LEVEL > 1) cat('\njcu=', i, 'cc=', jcu.cc[i], 'stage3.pop=',pop[3,i], 'floater.pop=',pop[4,i],
            'no over cc=', max(pop[3,i] - jcu.cc[i],0) , '\n')


        # ------------------------------------------------------------
        # Do floater dispersal into attached JCU
        # ------------------------------------------------------------

        # first check if there any floaters associated with the JCU and if any
        # territories have become available (i.e. due to mortality the
        # population is now below K)
        num.floaters.moving.into.jcu <- determine.floaters.moving.to.jcu(jcu.cc[source.jcu], 
                                            pop['stage3', source.jcu], pop['floaters', source.jcu])

        if(DEBUG.LEVEL > 1 & num.floaters.moving.into.jcu >0 ) {
            cat( '\n***************', num.floaters.moving.into.jcu, 'floater(s) moving into their attached JCU\n') 
            #browser()
        }

        # update the population matrix
        pop['stage3', source.jcu] <- pop['stage3', source.jcu] + num.floaters.moving.into.jcu
        pop['floaters', source.jcu] <- pop['floaters', source.jcu] - num.floaters.moving.into.jcu
        
        
        # calculate the number of adult individuals above carrying capactiy, as they
        # are the ones we assume will disperse
		num.stage3.to.disperse <- max(pop['stage3',source.jcu] - jcu.cc[source.jcu], 0)
        num.floaters.to.disperse <- pop['floaters', source.jcu]
        total.to.disperse <- num.stage3.to.disperse + num.floaters.to.disperse


        # reduce the pop of the source jcu by the number that disperse
        pop['stage3',source.jcu] <- pop['stage3',source.jcu] - num.stage3.to.disperse
        pop['floaters',source.jcu] <- 0

		# if there are more stage 3 adults than the cc disperse them
		if( total.to.disperse > 0) {

            # ------------------------------------------------------------
            # Work out where they disperse to and how many die  on the way
            # ------------------------------------------------------------

            # choose the jcus that each jaguar will try and disperse to;
            # select each one randomly for now NOTE: here is where to change
            # things if we want to try something other than random choices of
            # where to disperse to

            # TODO: need to limit which JCUs they can reach based on
            # assumption of max dispersal distance
			dest.jcus <- sample(jcu.vec[-source.jcu], total.to.disperse, replace=TRUE )
			
			# determine which ones survive the dispersal
			survival.vec <- rbinom(n=total.to.disperse, size=1, prob=(1-disp.mort.mat[source.jcu, dest.jcus]) )
			dest.jcus.surviving <- dest.jcus[which(survival.vec==1)]
			num.die.dispersing <- length(which(survival.vec==0))

			# determine the total number going to each jcu
			# Note: rle: Run Length Encoding to compute the lengths and values of runs of equal values 
			#            in a vector
			rle.cts <- rle(sort(dest.jcus.surviving))
			receive.cts <- rep(0, length(jcu.vec))
			receive.cts[rle.cts$values] <-  rle.cts$lengths

            # -------------------------------------------------------------
            # Of the individuals arriving in a patch, work out how many can
            # have territories (ie if the patch is not yet at it's carrying
            # capacity), and how many become floaters (assocuated with the
            # patch but not having a terriotiry)
            # -------------------------------------------------------------

            # for each JCU, determine how many (if any) home ranges are left to fill
            # before reaching K
            number.of.home.ranges.left <-  jcu.cc - pop['stage3',]
            number.of.home.ranges.left[number.of.home.ranges.left<=0] <- 0

            # determine how mant of the arriving individuals become floaters
            receive.cts.floaters <- receive.cts - number.of.home.ranges.left
            receive.cts.floaters[receive.cts.floaters<0] <- 0 

            # determine how many arriving individuals get territories (the
            # remainder that don't become floaters)
            receive.cts.stage3 <- receive.cts - receive.cts.floaters
            
			# update the population
            #pop[3,] <- pop['stage3',] + receive.cts
            pop['stage3',] <- pop['stage3',] + receive.cts.stage3
			pop['floaters',] <- pop['floaters',] + receive.cts.floaters

			# track some dispersal stats
			dispersal.ctr <- dispersal.ctr + total.to.disperse
			dispersal.mort.ctr <- dispersal.mort.ctr + num.die.dispersing

			
		}


	}

	
	if(DEBUG.LEVEL>0) cat(' (Num disp:', dispersal.ctr, 'mort in disp:', dispersal.mort.ctr, ')')
	
	return (pop)
}

#---------------------------------------------------------------

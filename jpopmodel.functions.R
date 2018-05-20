
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
    		# pre.truncated.mean of 1.06 corresponds to truncated.mean of 1.62
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

apply.dispersal <- function(pop, jcu.cc, disp.mort.mat) {

	if(DEBUG.LEVEL > 1) {cat( '\n  In dispersal function\n'); show(pop) }

	jcu.vec <- 1:length(jcu.cc)
	dispersal.ctr <- 0
	dispersal.mort.ctr <- 0

    # vector to hold the total number arricing at each patch after all
    # dispersal is done
    receive.cts <- rep(0, length(jcu.vec))
	for( source.jcu in jcu.vec ){

        # print out current JCU
        i<- source.jcu
        if(DEBUG.LEVEL > 1) cat('\njcu=', i, 'cc=', jcu.cc[i], 'stage3.pop=',pop[3,i], 'floater.pop=',pop[4,i],
            'no over cc=', max(pop[3,i] - jcu.cc[i],0) , '\n')


        #------------------------------------------------------------
        # Do floater dispersal into attached JCU
        #------------------------------------------------------------

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
      
			receive.cts[rle.cts$values] <-  receive.cts[rle.cts$values] + rle.cts$lengths

            
            
			
		} else {
            num.die.dispersing <- 0
    }

	}

      

        # -------------------------------------------------------------
        # Of the individuals arriving in each JCU, work out how many can
        # have territories (ie if the JCU is not yet at it's carrying
        # capacity), and how many become floaters (assocuated with the
        # patch but not having a terriotiry)
        # -------------------------------------------------------------


    # for each JCU, determine how many (if any) home ranges are left to fill
    # # before reaching K
    number.of.home.ranges.left <-  jcu.cc - pop['stage3',]
    number.of.home.ranges.left[number.of.home.ranges.left<=0] <- 0

    # determine how many of the arriving individuals become floaters
    receive.cts.floaters <- receive.cts - number.of.home.ranges.left
    receive.cts.floaters[receive.cts.floaters<0] <- 0 

    # determine how many arriving individuals get territories (the
    # remainder that don't become floaters)
    receive.cts.stage3 <- receive.cts - receive.cts.floaters
    
    # update the population
    pop['stage3',] <- pop['stage3',] + receive.cts.stage3
    pop['floaters',] <- pop['floaters',] + receive.cts.floaters

    # track some dispersal stats
    dispersal.ctr <- dispersal.ctr + total.to.disperse
    dispersal.mort.ctr <- dispersal.mort.ctr + num.die.dispersing


	if(DEBUG.LEVEL>0) cat(' (Num disp:', dispersal.ctr, 'mort in disp:', dispersal.mort.ctr, ')')

	return (pop)
}

#---------------------------------------------------------------

# Function to sample model input params to provide multiple realizatons from each expert

get.parameters <- function(Elicitation, JCUs, Distances, reps) {
	#generate ensemble
	Ensemble <- expand.grid(1:reps,1:5)
	Ensemble <- cbind(1:nrow(Ensemble),Ensemble)
	names(Ensemble) <- c("ID","REP","EXPERT")

	#simulate values for non-JCU specific parameters

	#BIRTH RATE
	BirthRate <- apply(X=Ensemble,MARGIN=1,FUN=function(X,Ensemble){rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="BIRTH_LOW",X["EXPERT"] + 1],
        b=Elicitation[Elicitation[,"PARAM"]=="BIRTH_HIGH",X["EXPERT"] + 1],
        c=Elicitation[Elicitation[,"PARAM"]=="BIRTH_BEST",X["EXPERT"] + 1])},
    Ensemble=Ensemble)

	#LITTER SIZE
	Litter <- apply(X=Ensemble,MARGIN=1,FUN=function(X,Ensemble){
        rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="LITTER_LOW",X["EXPERT"] + 1],
            b=Elicitation[Elicitation[,"PARAM"]=="LITTER_HIGH",X["EXPERT"] + 1],
            c=Elicitation[Elicitation[,"PARAM"]=="LITTER_BEST",X["EXPERT"] + 1])},
        Ensemble=Ensemble)

	#DISPERSAL MORTALITY
	DispM <- apply(X=Ensemble,MARGIN=1,FUN=function(X,Ensemble){
        rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="DISPM_LOW",X["EXPERT"] + 1],
            b=Elicitation[Elicitation[,"PARAM"]=="DISPM_HIGH",X["EXPERT"] + 1],
            c=Elicitation[Elicitation[,"PARAM"]=="DISPM_BEST",X["EXPERT"] + 1])},
        Ensemble=Ensemble)

	#MEAN DISPERSAL DISTANCE
	MeanD <- apply(X=Ensemble,MARGIN=1,FUN=function(X,Ensemble){rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MEAND_LOW",X["EXPERT"] + 1],
        b=Elicitation[Elicitation[,"PARAM"]=="MEAND_HIGH",X["EXPERT"] + 1],
        c=Elicitation[Elicitation[,"PARAM"]=="MEAND_BEST",X["EXPERT"] + 1])},
    Ensemble=Ensemble)

	#MAXIMUM DISPERSAL DISTANCE
	MaxD <- apply(X=Ensemble,MARGIN=1,FUN=function(X,Ensemble){rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MAXD_LOW",X["EXPERT"] + 1],
        b=Elicitation[Elicitation[,"PARAM"]=="MAXD_HIGH",X["EXPERT"] + 1],
        c=Elicitation[Elicitation[,"PARAM"]=="MAXD_BEST",X["EXPERT"] + 1])},
    Ensemble=Ensemble)

	#create non-JCU specific
	JCUInd_Params <- as.data.frame(cbind(BirthRate, Litter, DispM, MeanD, MaxD))
	names(JCUInd_Params) <- c("BirthR","LitS","DispM","MeanD","MaxD")

	#simulate values for JCU specific parameters

	#CARRYING CAPACITY
	K <- apply(X=Ensemble,MARGIN=1,FUN=function(X,Ensemble,JCUs){
        Size <- JCUs[,"Km2"];
        HR <- rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="HR_LOW",X["EXPERT"] + 1],
            b=Elicitation[Elicitation[,"PARAM"]=="HR_HIGH",X["EXPERT"] + 1],
            c=Elicitation[Elicitation[,"PARAM"]=="HR_BEST",X["EXPERT"] + 1]);
        return((Size*100)/HR)},
        Ensemble=Ensemble,JCUs=JCUs)

	#MORTALITY YOUNG
	MORTY <- apply(X=Ensemble,MARGIN=1,FUN=function(X,Ensemble,JCUs){
        Hunting <- JCUs[,"HunPress"];
        Mort <- ifelse(Hunting==1,rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MORTY1_LOW",X["EXPERT"] + 1],
            b=Elicitation[Elicitation[,"PARAM"]=="MORTY1_HIGH",X["EXPERT"] + 1],
            c=Elicitation[Elicitation[,"PARAM"]=="MORTY1_BEST",X["EXPERT"] + 1]),ifelse(Hunting==2,rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MORTY2_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="MORTY2_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="MORTY2_BEST",X["EXPERT"] + 1]),rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MORTY3_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="MORTY3_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="MORTY3_BEST",X["EXPERT"] + 1])));return(Mort)},Ensemble=Ensemble,JCUs=JCUs)

	#MORTALITY ADOLESCENT
	MORTA <- apply(X=Ensemble, MARGIN=1, FUN=function(X,Ensemble,JCUs){Hunting <- JCUs[,"HunPress"];Mort <- ifelse(Hunting==1,rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MORTA1_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="MORTA1_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="MORTA1_BEST",X["EXPERT"] + 1]),ifelse(Hunting==2,rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MORTA2_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="MORTA2_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="MORTA2_BEST",X["EXPERT"] + 1]),rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MORTA3_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="MORTA3_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="MORTA3_BEST",X["EXPERT"] + 1])));return(Mort)},Ensemble=Ensemble,JCUs=JCUs)

	#MORTALITY BREEDING AGE (ADULTS)
	MORTB <- apply(X=Ensemble,MARGIN=1,FUN=function(X,Ensemble,JCUs){Hunting <- JCUs[,"HunPress"];Mort <- ifelse(Hunting==1,rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MORTB1_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="MORTB1_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="MORTB1_BEST",X["EXPERT"] + 1]),ifelse(Hunting==2,rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MORTB2_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="MORTB2_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="MORTB2_BEST",X["EXPERT"] + 1]),rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MORTB3_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="MORTB3_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="MORTB3_BEST",X["EXPERT"] + 1])));return(Mort)},Ensemble=Ensemble,JCUs=JCUs)

	#INITIAL ADULT POPULATION SIZE
	START_POP <- JCUs[,"StartPop"]

	Output <- list(Ensemble, JCUInd_Params, K, MORTY, MORTA, MORTB, START_POP)

	names(Output) <- c("ENSEMBLE.LIST","FIXED.PARAMS","K","MORT.STAGE1","MORT.STAGE2","MORT.STAGE3","INITIAL.POP.STAGE3")

	return(Output)

}

#---------------------------------------------------------------

#function to generate the parameter realisations from the elicited data
#at the moment this usues a triangle distribution between the min and max
#with the best guess as the peak

get.parameters <- function(Elicitation,JCUs,Distances,Reps) {
	#generate ensemble of experts
	Ensemble <- expand.grid(1:Reps,1:5)
	Ensemble <- cbind(1:nrow(Ensemble),Ensemble)
	names(Ensemble) <- c("ID","REP","EXPERT")

	#simulate values for non-JCU specific parameters

	#BIRTH RATE
	BirthRate <- apply(X=Ensemble,MARGIN=1,FUN=function(X,Ensemble){rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="BIRTH_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="BIRTH_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="BIRTH_BEST",X["EXPERT"] + 1])},Ensemble=Ensemble)

	#LITTER SIZE
	Litter <- apply(X=Ensemble,MARGIN=1,FUN=function(X,Ensemble){rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="LITTER_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="LITTER_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="LITTER_BEST",X["EXPERT"] + 1])},Ensemble=Ensemble)

	#DISPERSAL MORTALITY
	DispM <- apply(X=Ensemble,MARGIN=1,FUN=function(X,Ensemble){rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="DISPM_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="DISPM_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="DISPM_BEST",X["EXPERT"] + 1])},Ensemble=Ensemble)

	#MEAN DISPERSAL DISTANCE
	MeanD <- apply(X=Ensemble,MARGIN=1,FUN=function(X,Ensemble){rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MEAND_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="MEAND_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="MEAND_BEST",X["EXPERT"] + 1])},Ensemble=Ensemble)

	#MAXIMUM DISPERSAL DISTANCE
	MaxD <- apply(X=Ensemble,MARGIN=1,FUN=function(X,Ensemble){rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MAXD_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="MAXD_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="MAXD_BEST",X["EXPERT"] + 1])},Ensemble=Ensemble)

	#create non-JCU specific parameters
	JCUInd_Params <- as.data.frame(cbind(BirthRate,Litter,DispM,MeanD,MaxD))
	names(JCUInd_Params) <- c("BirthR","LitS","DispM","MeanD","MaxD")

	#simulate values for JCU specific parameters

	#CARRYING CAPACITY
	K <- apply(X=Ensemble,MARGIN=1,FUN=function(X,Ensemble,JCUs){Size <- JCUs[,"Km2"];HR <- rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="HR_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="HR_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="HR_BEST",X["EXPERT"] + 1]) ;return((Size*100)/HR)},Ensemble=Ensemble,JCUs=JCUs)

	#MORTALITY YOUNG
	MORTY <- apply(X=Ensemble,MARGIN=1,FUN=function(X,Ensemble,JCUs){Hunting <- JCUs[,"HunPress"];Mort <- ifelse(Hunting==1,rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MORTY1_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="MORTY1_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="MORTY1_BEST",X["EXPERT"] + 1]),ifelse(Hunting==2,rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MORTY2_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="MORTY2_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="MORTY2_BEST",X["EXPERT"] + 1]),rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MORTY3_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="MORTY3_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="MORTY3_BEST",X["EXPERT"] + 1])));return(Mort)},Ensemble=Ensemble,JCUs=JCUs)

	#MORTALITY ADOLESCENT
	MORTA <- apply(X=Ensemble,MARGIN=1,FUN=function(X,Ensemble,JCUs){Hunting <- JCUs[,"HunPress"];Mort <- ifelse(Hunting==1,rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MORTA1_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="MORTA1_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="MORTA1_BEST",X["EXPERT"] + 1]),ifelse(Hunting==2,rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MORTA2_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="MORTA2_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="MORTA2_BEST",X["EXPERT"] + 1]),rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MORTA3_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="MORTA3_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="MORTA3_BEST",X["EXPERT"] + 1])));return(Mort)},Ensemble=Ensemble,JCUs=JCUs)

	#MORTALITY BREEDING AGE (ADULTS)
	MORTB <- apply(X=Ensemble,MARGIN=1,FUN=function(X,Ensemble,JCUs){Hunting <- JCUs[,"HunPress"];Mort <- ifelse(Hunting==1,rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MORTB1_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="MORTB1_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="MORTB1_BEST",X["EXPERT"] + 1]),ifelse(Hunting==2,rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MORTB2_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="MORTB2_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="MORTB2_BEST",X["EXPERT"] + 1]),rtriangle(n=1,a=Elicitation[Elicitation[,"PARAM"]=="MORTB3_LOW",X["EXPERT"] + 1],b=Elicitation[Elicitation[,"PARAM"]=="MORTB3_HIGH",X["EXPERT"] + 1],c=Elicitation[Elicitation[,"PARAM"]=="MORTB3_BEST",X["EXPERT"] + 1])));return(Mort)},Ensemble=Ensemble,JCUs=JCUs)

	#INITIAL ADULT POPULATION SIZE
	START_POP <- JCUs[,"StartPop"]

	Output <- list(Ensemble,JCUInd_Params,t(K),t(MORTY),t(MORTA),t(MORTB),START_POP,Distances)
	names(Output) <- c("ENSEMBLE.LIST","FIXED.PARAMS","K","MORT.STAGE1","MORT.STAGE2","MORT.STAGE3","INITIAL.POP.STAGE3","DISTANCES")

	return(Output)
}

#---------------------------------------------------------------

#function that takes the output from get.parameters() and formats the parameters
#for calling run.jpop.model using apply on the $ENSEMBLE.LIST matrix from the
#get.parameters() output
#Ensemble = ensemble, i.e. $ENSEMBLE.LIST from the output of get.parameters()
#Params.List = output (a list) from get.parameters()
#Reps = number of stochastic replicates

run.pop.model.apply <- function(Ensemble, Params.List, Years, Reps) {

	#get expert information
	expert.ID <- Params.List$ENSEMBLE.LIST[Ensemble[1],"EXPERT"]
	expert.realization <- Params.List$ENSEMBLE.LIST[Ensemble[1],"REP"]

	#create matrix to hold initial population size parameters
	initial.population <- matrix(0,nrow = 4, ncol = ncol(Params.List$K))
	colnames(initial.population) <- paste0( 'jcu', 1:ncol(Params.List$K)) # set column names
	rownames(initial.population) <- c(paste0( 'stage', 1:3), 'floaters' ) # set the row names

	#populate initial population size values
	#note that for now we assume that these are all adults
	initial.population[3,] <- Params.List$INITIAL.POP.STAGE3

	#create other carrying capacity, mortality, and birth rate parameters
	jcu.attributes <- matrix(0, ncol = 6, nrow = ncol(Params.List$K))
	colnames(jcu.attributes) <- c("K", "mortality.stage1", "mortality.stage2", "mortality.stage3", "mortality.floaters", "birth.rate.mean") # set column names

	#populate matrix
	#get K
	jcu.attributes[,"K"] <- Params.List$K[Ensemble[1],]
	#mortality
	jcu.attributes[,"mortality.stage1"] <- Params.List$MORT.STAGE1[Ensemble[1],]
	jcu.attributes[,"mortality.stage2"] <- Params.List$MORT.STAGE2[Ensemble[1],]
	jcu.attributes[,"mortality.stage3"] <- Params.List$MORT.STAGE3[Ensemble[1],]
	jcu.attributes[,"mortality.floaters"] <- Params.List$MORT.STAGE3[Ensemble[1],]
	#birth rate
	jcu.attributes[,"birth.rate.mean"] <- Params.List$FIXED.PARAMS[Ensemble[1],"BirthR"]
	#create data frame
	jcu.attributes <- as.data.frame(jcu.attributes)

	#NOTE NEED TO DEAL WITH LITTER SIZE - NOT DEALT WITH THIS YET (JR - 20/5/18)

	#dispersal mortality
	#note that here for now we make some simple assumptions about the dispersal process
	#we assume that if the distance between JCUs is > the maximum dispersal distance
	#then dispersal mortality = 1, otherwise we assume that dispersal mortality is
	#= elicited dispersal mortality probability per 10 km ^ (distance_km / 10).

	#NEED TO CHECK WE ARE USING THE RIGHT DISPERSAL PARAMETERS - CURRENTLY USING ADULT PARAMETERS (JR 20/5/18)
	#NEED TO DISCUSS ASSUMPTIONS MADE HERE ABOUT DISPERSAL - ARE WE HAPPY WITH THIS? (JR 20/5/18)

	#calculate dispersal mortality
	disp.mort.matrix <- Params.List$FIXED.PARAMS[Ensemble[1],"DispM"] ^ (Params.List$DISTANCES / 10)
	#account for maximum dispersal distance - set mortality to 1 if greater than maximum dispersal
	disp.mort.matrix <- ifelse(Params.List$DISTANCES > Params.List$FIXED.PARAMS[Ensemble[1],"MaxD"],1,disp.mort.matrix)

	#run population model
	run.jpop.model(expert.ID = expert.ID, expert.realization = expert.realization, num.stoch.realizatons = Reps, initial.pop = initial.population, jcu.att = jcu.attributes, disp.mort.mat = disp.mort.matrix, num.time.steps = Years)

}

#---------------------------------------------------------------

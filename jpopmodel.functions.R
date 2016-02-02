
#---------------------------------------------------------------

age.population <- function( cur.pop, num.stages) {

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

reproduce <- function( cur.pop, stages, birth.rate ) {

    tmp <- cur.pop
    s1 <- cur.pop[1,]
    pop <- rep(0, length(s1) )

    for( i in 1:length(s1) ) {
        pop[i] <- s1[i] + rbinom(n=1, size=cur.pop[stages,i], prob=birth.rate[i])
    }
    
    cur.pop[1,] <- pop
    

    return ( cur.pop  )
}

#---------------------------------------------------------------

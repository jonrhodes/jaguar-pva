library("triangle")

source('jpopmodel.functions.R')


Elicitation <- read.csv("input_data/elicited_parameters.csv",header=T)
JCUs <- read.csv("input_data/jcu_hunting_start_pop.csv",header=T)
Distance <- read.csv("input_data/distmatrix.csv",header=F)

Test <- get.parameters(Elicitation,JCUs,Distances,reps=5)

#TestFun <- function(X,Extract1,Extract2){
#preprocessing then call population model

#print birth rates to screen
# show(Extract1[X["ID"],"BirthR"])

# show(Extract2[,X["ID"]])

# pop_model(.....)

# return("success")

# }

# apply(X=Test$ENSEMBLE.LIST,MARGIN=1,FUN=TestFun,Extract1=Test$FIXED.PARAMS,Extract2=Test$K)


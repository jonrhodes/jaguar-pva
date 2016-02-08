# jaguar-pva
A PVA model for Barzillian Jagual in the Atlantic forest region.

To run:

1. Download the code from https://github.com/ascelin/jaguar-pva
2. Start R and make the working directory, the directory that contains
   the jaguar-pva code
3. To run the mode execute  source( 'jpopmodel.R') in R
4. To analyses the model, run results run source( 'analyse.results.R' ) in R


Current state of the code 

* Currently has no dispersal (this will be the next thing to be added)

* There are 3 life stages and only the 3rd stage can reproduce

* Assuming a single mortality facor for each JCU applies to each life
  stage

* For now reproduction is simply modelled follows: each individual in
  stage 3 has a predefined propability of having 1 offspring. This is
  all contained in the "reproduce" function.
  
Assumptions of the model

* Female only

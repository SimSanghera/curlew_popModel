######  Simple Curlew Population Viability Analysis #####

# Model adapted from Shoemaker and Oppel
# This scenario is to see if the model can match actual trends based of 
# curlew decline from 1987 - 2013.
# Using the estimated pop size from 1987 as the initial population size
# Predicting on a year-by-year basis, we will see how close the current
# setup is. 
# NOTE: literature has already stated that specific events caused shifts
# in survival and productivity rates. The current pop model utilises a
# single prod rate and survival rates, with only a minimal amount of 
# annual stochasticity built in - based on random catastrophes and yearly
# breeding failures.
# Whilst it would be possible to change the productivity rates for 
# certain years, this may not be helpful for the model in predicting future
# populations as we currently don't know the impact of any future rule changes
# or environmental events.


#----- Load Libraries -----
libs <- c("popbio",
          "doParallel",
          "foreach",
          "tidyverse",
          "markdown",
          "rmarkdown",
          "knitr",
          "data.table")

lapply(libs,
       require,
       character.only = TRUE)


#-----  Simulate Population Trajectory Across a Range of Parameters -----

# Specify parameters for the model
# Population estimates are counts of breeding pairs - thus assumes
# 0 chicks and juveniles at start of sequence
# NI.pop.sizes    <- seq(250, 560, 10)      # NI pop range from 250 - 500 
# # with 10 pair intervals
# 
# ROI.pop.sizes   <- seq(100, 200, 10)      # ROI pop range
# 
# I.popsizes      <- seq(100, 560, 20)      # all ireland pop size

# Pop Sizes
pop.1987 <- seq(3800, 6250, 200)
pop.1999 <- seq(1240, 2928, 200)
pop.2013 <- seq(252, 783, 200)

S.ad            <- seq(0.60, 0.92, 0.02)  # adult survival
# rates obtained from literature
# using UK wide rates
# NEED TO OBTAIN EU RATES

S.jw            <- seq(0.14, 0.49, 0.02)  # Juvenile survival - chick to
# fledging. UK wide. GET EU Numbers

S.jc            <- seq(0.20, 0.55, 0.02)  # SO suggested fledge survival rates
# (prod) were too low, so try with 
# a separate age class of released 
# chicks, with assumed higher
# survival

prod.rates     <- seq(0.14, 0.50, 0.02)  # productivity = number of fledged
# chicks per pair

released        <- seq(10, 50, 10)        # the number of birds to be
# introduced into the pop each year

intro.years     <- seq(1, 5, 1)            # the number of years in which
# birds will be introduced


# Specify parameters for Pop Viability Analysis
n.years   <- 30                           # number of years to simulate
n.reps    <- 500                          # number of simulations
K         <- 20000                        # carrying capacity of site
capt.fail <- 0.10                         # proportion of captive birds failed
# this is prior to release in the
# wild


# Specify parameters for stochastic components of the PVA
SD.lambda   <- 0.10                       # estimate lambda for pop growth rate
                                          # SD allows for building 
                                          # uncertainty into the model

phil.p      <- 0.05                       # probability that birds don't return
                                          # to release/breeding site 
                                          # currently unknown, so make estimate

return.p    <- 0.95                       # this might not be needed as it is 
                                          # a philopatric trend, so once a 
                                          # bird has bred at a site, it shows
                                          # fidelity

catastrophe.p <- 0.001                    # simulate potential for collapse of 
                                          # population due to catastrophic 
                                          # event
catastrophe.i <- 0.25                     # impact of catastrophe on pop
                                          # i.e. prop of surviving birds
                                          # should we differentiate impact 
                                          # between adults and juvs?

breed.fail    <- 0.1                      # prob the breeding season is a 
                                          # complete failure
                                          # it might be possible to find this
                                          # in lit and assign real number

extinct.cutoff  <- c(1, 2, 5, 10)         # below this number, pop is extinct
acceptable.risk <- 0.1                    # risk considered acceptable to 
                                          # managers

sex.ratio <- c(1)                       # in first model assume 50:50 but 
                                        # could be represented by a range


#-----  Defining Functions for Stochastic PVA -----

# Defining the structure of the population model with adult & juvenile stages
# Need to add other age stages representing survival
#   - egg to hatch
#   - 1 year breeders
#   - S.jw (and S.jr) is the same as prod.rate as we are assuming same values
# And prod.rate is the proportion of chicks that make it to next stage

bird.matrix <- expression(
  0,    prod.rates*0.5,
  S.jw*phil.p, S.ad
)



# Density-dependent Ricker model to simulate population growth

Ricker <- function(prev.abund) {
  
  # function for calculating next year abundance
  # includes environmental stochasticity
  
  prev.abund * exp(log(rnorm(1,                 # multiply prev.abund by 
                             R.max,             # the exp of the log of 
                             SD.lambda))        # a single value of growth rate
                   * (1 - (prev.abund/K)))      # from a normal distribution
  # using the SD of lambda
}


R.max <- 1            # set max pop growt rate


#-----  Population Projection Function  -----
PVA <- function(
  # arguments required 
  n.reps,
  n.years,
  init.N,
  K,
  catastrophe.p,
  catastrophe.i,
  F.loop,
  S.jw.loop,
  S.ad.loop,
  philopatry,
  breed.fail
) {
  
  # Create array to store population values
  # rows = years, cols = reps
  Pop.Array <- array(0,
                     dim = c((n.years + 1),
                             n.reps))
  
  # Start loop of reps
  for (rep in 1:n.reps) {
    
    # Initial abundance minus number of head-started birds
    #   but we will have a sequence of head-starting scenarios
    #   so won't need to select from a poisson distribution
    Pop.Array[1, rep]   <- init.N * 2 
    
    
    # Loop through the years
    for (y in 2:(n.years + 1)) {
      
      # create matrix with vital rates
      # Ruined breeding season due to whatever...
      # Possibly count number of poor breeding years in survey data and 
      #   obtain more accurate probability for this occurring
      # F.loop comes from simul.in table = fecundity
      F.Year <- ifelse(rbinom(1, 1, breed.fail) == 1,
                       F.loop * 0.5,
                       F.loop)
      
      # Create vital rates
      bird.vr <- list(prod.rates  = F.Year,
                      S.ad    = S.ad.loop,
                      S.jw    = S.jw.loop,
                      phil.p  = philopatry)
      
      
      # Create Leslie Matrix
      a <- matrix(sapply(bird.matrix,
                         eval,
                         bird.vr,
                         NULL),
                  nrow = sqrt(length(bird.matrix)),
                  byrow = TRUE)
      
      # Starting population vector
      # 0 chicks in our model as all literature reports breeding pairs
      pop.size  <- c(0,
                     (init.N * 2))  # Multiplied by 2 to 
                                    #   obtain individuals
      
      # Get projections for years
      projections <- pop.projection(a,
                                    n = pop.size,
                                    iterations = 30)
      
      # Stochastic Ricker Model
      R.max <- projections$lambda     # growth rate calculated by model
      
      # Calculate abundance
      next.year <- max(0,
                       trunc(Ricker(Pop.Array[y - 1, rep])))
      
      # Catastrophe
      if(runif(1) < catastrophe.p) next.year <- next.year * catastrophe.i 
      
      # Catastrophe with head-starting
      # if(y < intro.years) next.year <- next.year + init.N
      
      # Put new abundance in array
      Pop.Array[y, rep] <- next.year
      
    }   # end year
    
  }     # end rep
  
  return(Pop.Array)
  
}


#---- test  ----
testPva <- PVA(
  # arguments required 
  n.reps,
  n.years,
  init.N <- simul.in[1, 1],
  K,
  catastrophe.p,
  catastrophe.i,
  F.loop <- simul.in[1, 4],
  S.jw.loop <- simul.in[1, 3],
  S.ad.loop <- simul.in[1, 2],
  philopatry <- simul.in[1, 5],
  breed.fail
) 

#-----  Calculate Proportion of Sims where Species Go Extinct -----
extinction.bysim  <- function(simdata,
                              threshold)  {
  sum(apply(default,
            2,
            function(t) min(t) < threshold)) / ncol(simdata)
  # extinction is defined as < threshold number of birds
}


#-----  Loop Over Combinations of Demography  -----
# Table of all combinations of demographic parameters

simul.in <- expand.grid(init.N    = pop.1987,
                        S.ad.loop = S.ad,
                        S.jw.loop = S.jw,
                        F.loop    = prod.rates,
                        philopatry = phil.p) # for release, add intro years

dim(simul.in)
names(simul.in) <- c("pop.size",
                     "S.ad",
                     "S.jw",
                     "Productivity Rate (F*2)")


# Test values
# 12747
# F.loop      <- simul.in[12747, 4]
# S.ad.loop   <- simul.in[12747, 2]
# S.jw.loop   <- simul.in[12747, 3]
# 
# 
# #https://dcl-prog.stanford.edu/purrr-basics.html
# #https://dcl-prog.stanford.edu/purrr-parallel.html
# #https://www.r-bloggers.com/2021/09/running-r-code-for-all-combinations-of-some-parameters-with-lapply-karate/
# 
# init.N <- 2000
# 
# # Create small combo grid
# Prod.rate <- seq(0.2, 0.5, 0.1)
# S.ad.loop <- seq(0.8, 0.9, 0.2)
# S.jw.loop <- seq(0.2, 0.6, 0.2)
# 
# simul.in <- expand.grid(F.loop = Prod.rate,
#                         S.ad.loop = S.ad.loop,
#                         S.jw.loop = S.jw.loop)

i <- 1:nrow(simul.in)

system.time(
  sim.out <- with(simul.in,
                lapply(i,
                       function(j){PVA(n.reps,
                                       n.years,
                                       init.N[j],
                                       K,
                                       catastrophe.p,
                                       catastrophe.i,
                                       F.loop[j],
                                       S.jw.loop[j],
                                       S.ad.loop[j],
                                       philopatry[j],
                                       breed.fail)}))
)



#-----  Try Parallel Processing -----

# Boost speed using parallel processing 
cl <- makeCluster(8)
registerDoParallel(cl, cores = 8)


# Start loop
system.time(
  
  sim.out <- foreach(s = c(1:dim(simul.in)[1]),
                     .packages = "popbio") %dopar% {
        # change to orginal code - sim.out is now a list, every rep is
        #   its own entry. removed .combine = cbind from above
                       
                       # start pva run
                       init.N <- simul.in[s, 1]
                         default <- PVA(n.reps,
                                        n.years,
                                        init.N,
                                        K, 
                                        catastrophe.p,
                                        catastrophe.i,
                                        F.loop      = simul.in[s, 4],
                                        S.jw.loop   = simul.in[s, 3],
                                        S.ad.loop   = simul.in[s, 2],
                                        philopatry  = simul.in[s, 5],
                                        breed.fail)
                       
                     }
  
)

# Stop parallel processing
stopCluster(cl)

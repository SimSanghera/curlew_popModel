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

# Pop Sizes NI
pop.1987 <- seq(3800, 6250, 200)
pop.1999 <- seq(1240, 2928, 200)
pop.2013 <- seq(252, 783, 200)

# Pop size ROI
ROI.pop.1985 <- seq(5000, 9000, 500)


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
n.reps    <- 500                         # number of simulations
K         <- 20000                        # carrying capacity of site
ROI.K     <- 25000
capt.fail <- 0.10                         # proportion of captive birds failed
# this is prior to release in the
# wild


# Specify parameters for stochastic components of the PVA
SD.lambda   <- 0.10                       # estimate lambda for pop growth rate
                                          # SD allows for building 
                                          # uncertainty into the model

phil.p      <- seq(0.30, 0.90, 0.5)       # probability that birds don't return
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
      # v3 - using newly calculated chick and adult numbers
      if (y == 2){
        pop.size <- c(0, (init.N*2))
      } else {
        pop.size <- c(projections$stage.vectors[, y-1])
      }
     
      # pop.size  <- c(0,
      #                (init.N * 2))  # Multiplied by 2 to 
      #                               #   obtain individuals
      
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


prev.abund <- Pop.Array[y-1, rep]
newpop <- exp(log(rnorm(1, R.max, SD.lambda))*(1-(prev.abund/K)))
newabund <- prev.abund*newpop
next.year <- max(0, trunc(newabund))
Pop.Array[y, rep] <- next.year

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


write.table(simul.in, "data/base_mod_1987_present_simulatons.csv",
            sep = ",",
            row.names = FALSE)
NI.sim.params <- read.csv("data/base_mod_1987_present_simulatons.csv",
                          header =  TRUE)

# Test values
# 12747
# init.N      <- simul.in[12747, 1]
# F.loop      <- simul.in[12747, 4]
# S.ad.loop   <- simul.in[12747, 2]
# S.jw.loop   <- simul.in[12747, 3]
# philopatry  <- simul.in[12747, 5]
# 
# init.N      <- simul.in[1, 1]
# F.loop      <- simul.in[1, 4]
# S.ad.loop   <- simul.in[1, 2]
# S.jw.loop   <- simul.in[1, 3]
# philopatry  <- simul.in[1, 5]
# 
# init.N      <- simul.in[27, 1]
# F.loop      <- simul.in[27, 4]
# S.ad.loop   <- simul.in[27, 2]
# S.jw.loop   <- simul.in[27, 3]
# philopatry  <- simul.in[27, 5]
# 
# init.N      <- simul.in[53, 1]
# F.loop      <- simul.in[53, 4]
# S.ad.loop   <- simul.in[53, 2]
# S.jw.loop   <- simul.in[53, 3]
# philopatry  <- simul.in[53, 5]

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
                       #init.N <- simul.in[s, 1]
                       # shift into function argument
                       default <- PVA(n.reps,
                                      n.years,
                                      init.N      = simul.in[s, 1],
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

names.NI.sim.out <- rep(as.character(c(pop.1987)), (151164/13))




# Write to save file
saveRDS(sim.out,
        "data/base_mod_1987_present_sims_3.RDS")
#ver 3 = breed.fail 0,4

#-----  Plotting  -----
NI.sim.out <- readRDS("data/base_mod_1987_present_sims_3.RDS")

names(NI.sim.out) <- names.NI.sim.out

years <- seq(1987, 2017, 1)

NI.sim.out.3800 <- NI.sim.out[names(NI.sim.out) == "3800"]
NI.sim.out.3800 <- do.call(cbind, NI.sim.out.3800)
NI.sim.out.3800.mean <- rowMeans(NI.sim.out.3800)
NI.sim.out.3800.lcl  <- apply(NI.sim.out.3800, 1, quantile, probs = 0.1)
NI.sim.out.3800.ucl  <- apply(NI.sim.out.3800, 1, quantile, probs = 0.9)
which(NI.sim.out.3800 == min(NI.sim.out.3800),
      arr.ind = TRUE)

which(NI.sim.out.3800 == max(NI.sim.out.3800),
      arr.ind = TRUE)

bettermatch <- which(NI.sim.out.3800[31 ,] < 500,
      arr.ind = TRUE)

NI.sim.out.3800.match <- NI.sim.out.3800[, bettermatch]
NI.sim.out.3800.match.mean <- rowMeans(NI.sim.out.3800.match)
NI.sim.out.3800.match.lcl <- apply(NI.sim.out.3800.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.3800.match.ucl <- apply(NI.sim.out.3800.match, 1,
                                   quantile, probs = 0.9)

NI.sim.out.4000 <- NI.sim.out[names(NI.sim.out) == "4000"]
NI.sim.out.4000 <- do.call(cbind, NI.sim.out.4000)
NI.sim.out.4000.mean <- rowMeans(NI.sim.out.4000)
NI.sim.out.4000.lcl  <- apply(NI.sim.out.4000, 1, quantile, probs = 0.1)
NI.sim.out.4000.ucl  <- apply(NI.sim.out.4000, 1, quantile, probs = 0.9)

NI.sim.out.4200 <- NI.sim.out[names(NI.sim.out) == "4200"]
NI.sim.out.4200 <- do.call(cbind, NI.sim.out.4200)
NI.sim.out.4200.mean <- rowMeans(NI.sim.out.4200)
NI.sim.out.4200.lcl  <- apply(NI.sim.out.4200, 1, quantile, probs = 0.1)
NI.sim.out.4200.ucl  <- apply(NI.sim.out.4200, 1, quantile, probs = 0.9)

NI.sim.out.4400 <- NI.sim.out[names(NI.sim.out) == "4400"]
NI.sim.out.4400 <- do.call(cbind, NI.sim.out.4400)
NI.sim.out.4400.mean <- rowMeans(NI.sim.out.4400)
NI.sim.out.4400.lcl  <- apply(NI.sim.out.4400, 1, quantile, probs = 0.1)
NI.sim.out.4400.ucl  <- apply(NI.sim.out.4400, 1, quantile, probs = 0.9)

NI.sim.out.4600 <- NI.sim.out[names(NI.sim.out) == "4600"]
NI.sim.out.4600 <- do.call(cbind, NI.sim.out.4600)
NI.sim.out.4600.mean <- rowMeans(NI.sim.out.4600)
NI.sim.out.4600.lcl  <- apply(NI.sim.out.4600, 1, quantile, probs = 0.1)
NI.sim.out.4600.ucl  <- apply(NI.sim.out.4600, 1, quantile, probs = 0.9)

NI.sim.out.4800 <- NI.sim.out[names(NI.sim.out) == "4800"]
NI.sim.out.4800 <- do.call(cbind, NI.sim.out.4800)
NI.sim.out.4800.mean <- rowMeans(NI.sim.out.4800)
NI.sim.out.4800.lcl  <- apply(NI.sim.out.4800, 1, quantile, probs = 0.1)
NI.sim.out.4800.ucl  <- apply(NI.sim.out.4800, 1, quantile, probs = 0.9)

NI.sim.out.5000 <- NI.sim.out[names(NI.sim.out) == "5000"]
NI.sim.out.5000 <- do.call(cbind, NI.sim.out.5000)
NI.sim.out.5000.mean <- rowMeans(NI.sim.out.5000)
NI.sim.out.5000.lcl  <- apply(NI.sim.out.5000, 1, quantile, probs = 0.1)
NI.sim.out.5000.ucl  <- apply(NI.sim.out.5000, 1, quantile, probs = 0.9)

NI.sim.out.5200 <- NI.sim.out[names(NI.sim.out) == "5200"]
NI.sim.out.5200 <- do.call(cbind, NI.sim.out.5200)
NI.sim.out.5200.mean <- rowMeans(NI.sim.out.5200)
NI.sim.out.5200.lcl  <- apply(NI.sim.out.5200, 1, quantile, probs = 0.1)
NI.sim.out.5200.ucl  <- apply(NI.sim.out.5200, 1, quantile, probs = 0.9)

NI.sim.out.5400 <- NI.sim.out[names(NI.sim.out) == "5400"]
NI.sim.out.5400 <- do.call(cbind, NI.sim.out.5400)
NI.sim.out.5400.mean <- rowMeans(NI.sim.out.5400)
NI.sim.out.5400.lcl  <- apply(NI.sim.out.5400, 1, quantile, probs = 0.1)
NI.sim.out.5400.ucl  <- apply(NI.sim.out.5400, 1, quantile, probs = 0.9)

NI.sim.out.5600 <- NI.sim.out[names(NI.sim.out) == "5600"]
NI.sim.out.5600 <- do.call(cbind, NI.sim.out.5600)
NI.sim.out.5600.mean <- rowMeans(NI.sim.out.5600)
NI.sim.out.5600.lcl  <- apply(NI.sim.out.5600, 1, quantile, probs = 0.1)
NI.sim.out.5600.ucl  <- apply(NI.sim.out.5600, 1, quantile, probs = 0.9)

NI.sim.out.5800 <- NI.sim.out[names(NI.sim.out) == "5800"]
NI.sim.out.5800 <- do.call(cbind, NI.sim.out.5800)
NI.sim.out.5800.mean <- rowMeans(NI.sim.out.5800)
NI.sim.out.5800.lcl  <- apply(NI.sim.out.5800, 1, quantile, probs = 0.1)
NI.sim.out.5800.ucl  <- apply(NI.sim.out.5800, 1, quantile, probs = 0.9)

NI.sim.out.6000 <- NI.sim.out[names(NI.sim.out) == "6000"]
NI.sim.out.6000 <- do.call(cbind, NI.sim.out.6000)
NI.sim.out.6000.mean <- rowMeans(NI.sim.out.6000)
NI.sim.out.6000.lcl  <- apply(NI.sim.out.6000, 1, quantile, probs = 0.1)
NI.sim.out.6000.ucl  <- apply(NI.sim.out.6000, 1, quantile, probs = 0.9)

NI.sim.out.6200 <- NI.sim.out[names(NI.sim.out) == "6200"]
NI.sim.out.6200 <- do.call(cbind, NI.sim.out.6200)
NI.sim.out.6200.mean <- rowMeans(NI.sim.out.6200)
NI.sim.out.6200.lcl  <- apply(NI.sim.out.6200, 1, quantile, probs = 0.1)
NI.sim.out.6200.ucl  <- apply(NI.sim.out.6200, 1, quantile, probs = 0.9)


NI.sim.3800.df <- data.frame(years,
                         NI.sim.out.3800.mean,
                         NI.sim.out.3800.lcl,
                         NI.sim.out.3800.ucl)
names(NI.sim.3800.df) <- c("years", "mean", "lcl", "ucl")

# Combine into tibble

NI.sim.out.df <- tibble(years = as.factor(years), 
                         NI.sim.out.3800.mean,
                         NI.sim.out.3800.lcl,
                         NI.sim.out.3800.ucl,
                         NI.sim.out.4000.mean,
                         NI.sim.out.4000.lcl,
                         NI.sim.out.4200.ucl,
                         NI.sim.out.4400.mean,
                         NI.sim.out.4400.lcl,
                         NI.sim.out.4400.ucl,
                         NI.sim.out.4600.mean,
                         NI.sim.out.4600.lcl,
                         NI.sim.out.4600.ucl,
                         NI.sim.out.4800.mean,
                         NI.sim.out.4800.lcl,
                         NI.sim.out.4800.ucl,
                         NI.sim.out.5000.mean,
                         NI.sim.out.5000.lcl,
                         NI.sim.out.5000.ucl,
                         NI.sim.out.5200.mean,
                         NI.sim.out.5200.lcl,
                         NI.sim.out.5200.ucl,
                         NI.sim.out.5400.mean,
                         NI.sim.out.5400.lcl,
                         NI.sim.out.5400.ucl,
                         NI.sim.out.5600.mean,
                         NI.sim.out.5600.lcl,
                         NI.sim.out.5600.ucl,
                         NI.sim.out.5800.mean,
                         NI.sim.out.5800.lcl,
                         NI.sim.out.5800.ucl,
                         NI.sim.out.6000.mean,
                         NI.sim.out.6000.lcl,
                         NI.sim.out.6000.ucl,
                         NI.sim.out.6200.mean,
                         NI.sim.out.6200.lcl,
                         NI.sim.out.6200.ucl)

NI.init.N <- rep(pop, 3)

NI.sim.out.df.long <- data.frame(years = rep(years, times = 13),
                             initN = rep(pop.1987, each = 31),
                             mean = c(NI.sim.out.3800.mean,
                                          NI.sim.out.4000.mean,
                                          NI.sim.out.4200.mean,
                                          NI.sim.out.4400.mean,
                                          NI.sim.out.4600.mean,
                                          NI.sim.out.4800.mean,
                                          NI.sim.out.5000.mean,
                                          NI.sim.out.5200.mean,
                                          NI.sim.out.5400.mean,
                                          NI.sim.out.5600.mean,
                                          NI.sim.out.5800.mean,
                                          NI.sim.out.6000.mean,
                                          NI.sim.out.6200.mean),
                             lcl = c(NI.sim.out.3800.lcl,
                                         NI.sim.out.4000.lcl,
                                         NI.sim.out.4200.lcl,
                                         NI.sim.out.4400.lcl,
                                         NI.sim.out.4600.lcl,
                                         NI.sim.out.4800.lcl,
                                         NI.sim.out.5000.lcl,
                                         NI.sim.out.5200.lcl,
                                         NI.sim.out.5400.lcl,
                                         NI.sim.out.5600.lcl,
                                         NI.sim.out.5800.lcl,
                                         NI.sim.out.6000.lcl,
                                         NI.sim.out.6200.lcl),
                             ucl = c(NI.sim.out.3800.ucl,
                                         NI.sim.out.4000.ucl,
                                         NI.sim.out.4200.ucl,
                                         NI.sim.out.4400.ucl,
                                         NI.sim.out.4600.ucl,
                                         NI.sim.out.4800.ucl,
                                         NI.sim.out.5000.ucl,
                                         NI.sim.out.5200.ucl,
                                         NI.sim.out.5400.ucl,
                                         NI.sim.out.5600.ucl,
                                         NI.sim.out.5800.ucl,
                                         NI.sim.out.6000.ucl,
                                         NI.sim.out.6200.ucl)
                             )


write.table(NI.sim.out.df.long,
            "data/NI.sim.out.base.plot.1987.csv",
            sep = ",",
            row.names = FALSE)





NI.sim.out.plot <- ggplot(NI.sim.out.df.long) +
  facet_wrap(~initN) + 
  geom_line(aes(x = years,
                y = mean,
                colour = factor(initN)))

NI.sim.out.plot.ribbon <- NI.sim.out.plot + 
  geom_ribbon(aes(x = years,
                  ymin = lcl,
                  ymax = ucl),
              alpha = 0.3) + 
  geom_point(data = NI.pop.ests,
             aes(x = years,
                 y = popsize)) +
  xlab("Years") + 
  ylab("Curlew Population Size (individuals)") +
  ylim(0, 20000) + 
  theme(legend.position = "none",
        panel.background = element_rect(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())

NI.pop.ests <- data.frame(years = c(1987, 1999, 2013),
                          popsize = c((2*5000), (2*2091), (2*526)))


NI.sim.out.plot
NI.sim.out.plot.ribbon


NI.3800.plot <- ggplot(NI.sim.3800.df) + 
  geom_line(aes(x = years,
            y = mean))
NI.3800.ribbon <- NI.3800.plot + 
  geom_ribbon(aes(x = years,
              ymin = lcl,
              ymax = ucl),
              alpha = 0.1)



#----- Obtaining predictions that predict less than 500 curlew  -----
bettermatch <- which(NI.sim.out.3800[31 ,] < 500,
                     arr.ind = TRUE)
NI.sim.out.3800.match <- NI.sim.out.3800[, bettermatch]
NI.sim.out.3800.match.mean <- rowMeans(NI.sim.out.3800.match)
NI.sim.out.3800.match.lcl <- apply(NI.sim.out.3800.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.3800.match.ucl <- apply(NI.sim.out.3800.match, 1,
                                   quantile, probs = 0.9)

bettermatch <- which(NI.sim.out.4000[31 ,] < 500,
                     arr.ind = TRUE)
NI.sim.out.4000.match <- NI.sim.out.4000[, bettermatch]
NI.sim.out.4000.match.mean <- rowMeans(NI.sim.out.4000.match)
NI.sim.out.4000.match.lcl <- apply(NI.sim.out.4000.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.4000.match.ucl <- apply(NI.sim.out.4000.match, 1,
                                   quantile, probs = 0.9)

bettermatch <- which(NI.sim.out.4200[31 ,] < 500,
                     arr.ind = TRUE)
NI.sim.out.4200.match <- NI.sim.out.4200[, bettermatch]
NI.sim.out.4200.match.mean <- rowMeans(NI.sim.out.4200.match)
NI.sim.out.4200.match.lcl <- apply(NI.sim.out.4200.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.4200.match.ucl <- apply(NI.sim.out.4200.match, 1,
                                   quantile, probs = 0.9)

bettermatch <- which(NI.sim.out.4400[31 ,] < 500,
                     arr.ind = TRUE)
NI.sim.out.4400.match <- NI.sim.out.4400[, bettermatch]
NI.sim.out.4400.match.mean <- rowMeans(NI.sim.out.4400.match)
NI.sim.out.4400.match.lcl <- apply(NI.sim.out.4400.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.4400.match.ucl <- apply(NI.sim.out.4400.match, 1,
                                   quantile, probs = 0.9)

bettermatch <- which(NI.sim.out.4600[31 ,] < 500,
                     arr.ind = TRUE)
NI.sim.out.4600.match <- NI.sim.out.4600[, bettermatch]
NI.sim.out.4600.match.mean <- rowMeans(NI.sim.out.4600.match)
NI.sim.out.4600.match.lcl <- apply(NI.sim.out.4600.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.4600.match.ucl <- apply(NI.sim.out.4600.match, 1,
                                   quantile, probs = 0.9)

bettermatch <- which(NI.sim.out.4800[31 ,] < 500,
                     arr.ind = TRUE)
NI.sim.out.4800.match <- NI.sim.out.4800[, bettermatch]
NI.sim.out.4800.match.mean <- rowMeans(NI.sim.out.4800.match)
NI.sim.out.4800.match.lcl <- apply(NI.sim.out.4800.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.4800.match.ucl <- apply(NI.sim.out.4800.match, 1,
                                   quantile, probs = 0.9)

bettermatch <- which(NI.sim.out.5000[31 ,] < 500,
                     arr.ind = TRUE)
NI.sim.out.5000.match <- NI.sim.out.5000[, bettermatch]
NI.sim.out.5000.match.mean <- rowMeans(NI.sim.out.5000.match)
NI.sim.out.5000.match.lcl <- apply(NI.sim.out.5000.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.5000.match.ucl <- apply(NI.sim.out.5000.match, 1,
                                   quantile, probs = 0.9)

bettermatch <- which(NI.sim.out.5200[31 ,] < 500,
                     arr.ind = TRUE)
NI.sim.out.5200.match <- NI.sim.out.5200[, bettermatch]
NI.sim.out.5200.match.mean <- rowMeans(NI.sim.out.5200.match)
NI.sim.out.5200.match.lcl <- apply(NI.sim.out.5200.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.5200.match.ucl <- apply(NI.sim.out.5200.match, 1,
                                   quantile, probs = 0.9)

bettermatch <- which(NI.sim.out.5400[31 ,] < 500,
                     arr.ind = TRUE)
NI.sim.out.5400.match <- NI.sim.out.5400[, bettermatch]
NI.sim.out.5400.match.mean <- rowMeans(NI.sim.out.5400.match)
NI.sim.out.5400.match.lcl <- apply(NI.sim.out.5400.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.5400.match.ucl <- apply(NI.sim.out.5400.match, 1,
                                   quantile, probs = 0.9)

bettermatch <- which(NI.sim.out.5600[31 ,] < 500,
                     arr.ind = TRUE)
NI.sim.out.5600.match <- NI.sim.out.5600[, bettermatch]
NI.sim.out.5600.match.mean <- rowMeans(NI.sim.out.5600.match)
NI.sim.out.5600.match.lcl <- apply(NI.sim.out.5600.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.5600.match.ucl <- apply(NI.sim.out.5600.match, 1,
                                   quantile, probs = 0.9)

bettermatch <- which(NI.sim.out.5800[31 ,] < 500,
                     arr.ind = TRUE)
NI.sim.out.5800.match <- NI.sim.out.5800[, bettermatch]
NI.sim.out.5800.match.mean <- rowMeans(NI.sim.out.5800.match)
NI.sim.out.5800.match.lcl <- apply(NI.sim.out.5800.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.5800.match.ucl <- apply(NI.sim.out.5800.match, 1,
                                   quantile, probs = 0.9)

bettermatch <- which(NI.sim.out.6000[31 ,] < 500,
                     arr.ind = TRUE)
NI.sim.out.6000.match <- NI.sim.out.6000[, bettermatch]
NI.sim.out.6000.match.mean <- rowMeans(NI.sim.out.6000.match)
NI.sim.out.6000.match.lcl <- apply(NI.sim.out.6000.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.6000.match.ucl <- apply(NI.sim.out.6000.match, 1,
                                   quantile, probs = 0.9)

bettermatch <- which(NI.sim.out.6200[31 ,] < 500,
                     arr.ind = TRUE)
NI.sim.out.6200.match <- NI.sim.out.6200[, bettermatch]
NI.sim.out.6200.match.mean <- rowMeans(NI.sim.out.6200.match)
NI.sim.out.6200.match.lcl <- apply(NI.sim.out.6200.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.6200.match.ucl <- apply(NI.sim.out.6200.match, 1,
                                   quantile, probs = 0.9)


NI.sim.out.df.match <- data.frame(years = rep(years, times = 13),
                                 initN = rep(pop.1987, each = 31),
                                 mean = c(NI.sim.out.3800.match.mean,
                                          NI.sim.out.4000.match.mean,
                                          NI.sim.out.4200.match.mean,
                                          NI.sim.out.4400.match.mean,
                                          NI.sim.out.4600.match.mean,
                                          NI.sim.out.4800.match.mean,
                                          NI.sim.out.5000.match.mean,
                                          NI.sim.out.5200.match.mean,
                                          NI.sim.out.5400.match.mean,
                                          NI.sim.out.5600.match.mean,
                                          NI.sim.out.5800.match.mean,
                                          NI.sim.out.6000.match.mean,
                                          NI.sim.out.6200.match.mean),
                                 lcl = c(NI.sim.out.3800.match.lcl,
                                         NI.sim.out.4000.match.lcl,
                                         NI.sim.out.4200.match.lcl,
                                         NI.sim.out.4400.match.lcl,
                                         NI.sim.out.4600.match.lcl,
                                         NI.sim.out.4800.match.lcl,
                                         NI.sim.out.5000.match.lcl,
                                         NI.sim.out.5200.match.lcl,
                                         NI.sim.out.5400.match.lcl,
                                         NI.sim.out.5600.match.lcl,
                                         NI.sim.out.5800.match.lcl,
                                         NI.sim.out.6000.match.lcl,
                                         NI.sim.out.6200.match.lcl),
                                 ucl = c(NI.sim.out.3800.match.ucl,
                                         NI.sim.out.4000.match.ucl,
                                         NI.sim.out.4200.match.ucl,
                                         NI.sim.out.4400.match.ucl,
                                         NI.sim.out.4600.match.ucl,
                                         NI.sim.out.4800.match.ucl,
                                         NI.sim.out.5000.match.ucl,
                                         NI.sim.out.5200.match.ucl,
                                         NI.sim.out.5400.match.ucl,
                                         NI.sim.out.5600.match.ucl,
                                         NI.sim.out.5800.match.ucl,
                                         NI.sim.out.6000.match.ucl,
                                         NI.sim.out.6200.match.ucl)
)


NI.sim.out.match.plot <- ggplot(NI.sim.out.df.match) +
  facet_wrap(~initN) + 
  geom_line(aes(x = years,
                y = mean,
                colour = factor(initN)))

NI.sim.out.match.plot.ribbon <- NI.sim.out.match.plot + 
  geom_ribbon(aes(x = years,
                  ymin = lcl,
                  ymax = ucl),
              alpha = 0.3) + 
  geom_point(data = NI.pop.ests,
             aes(x = years,
                 y = popsize)) +
  xlab("Years") + 
  ylab("Curlew Population Size (individuals)") +
  ylim(0, 20000) + 
  theme(legend.position = "none",
        panel.background = element_rect(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())

NI.pop.ests <- data.frame(years = c(1987, 1999, 2013),
                          popsize = c((2*5000), (2*2091), (2*526)))



# Get combinations of parameters
bettermatch.3800 <- which(NI.sim.out.3800[31 ,] < 500,
                     arr.ind = TRUE)

NI.3800.match.params <- NI.sim.params[bettermatch.3800 ,]


# roi

ROI.Ricker <- function(prev.abund) {
  
  # function for calculating next year abundance
  # includes environmental stochasticity
  
  prev.abund * exp(log(rnorm(1,                 # multiply prev.abund by 
                             R.max,             # the exp of the log of 
                             SD.lambda))        # a single value of growth rate
                   * (1 - (prev.abund/ROI.K)))      # from a normal distribution
  # using the SD of lambda
}


#-----  Population Projection Function  -----
ROI.PVA <- function(
  # arguments required 
  n.reps,
  n.years,
  init.N,
  ROI.K,
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
      # v3 - using newly calculated chick and adult numbers
      if (y == 2){
        pop.size <- c(0, (init.N*2))
      } else {
        pop.size <- c(projections$stage.vectors[, y-1])
      }
      
      # pop.size  <- c(0,
      #                (init.N * 2))  # Multiplied by 2 to 
      #                               #   obtain individuals
      
      # Get projections for years
      projections <- pop.projection(a,
                                    n = pop.size,
                                    iterations = 30)
      
      # Stochastic Ricker Model
      R.max <- projections$lambda     # growth rate calculated by model
      
      # Calculate abundance
      next.year <- max(0,
                       trunc(ROI.Ricker(Pop.Array[y - 1, rep])))
      
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

# Combination of parameters
ROI.simul.in <- expand.grid(init.N    = ROI.pop.1985,
                        S.ad.loop = S.ad,
                        S.jw.loop = S.jw,
                        F.loop    = prod.rates,
                        philopatry = phil.p) # for release, add intro years

dim(simul.in)
names(simul.in) <- c("pop.size",
                     "S.ad",
                     "S.jw",
                     "Productivity Rate (F*2)")


write.table(ROI.simul.in, "data/ROI_base_mod_1985_present_simulatons.csv",
            sep = ",",
            row.names = FALSE)


#-----  Try Parallel Processing -----

breed.fail <- 0.5

# Boost speed using parallel processing 
cl <- makeCluster(8)
registerDoParallel(cl, cores = 8)


# Start loop
system.time(
  
  ROI.sim.out <- foreach(s = c(1:dim(ROI.simul.in)[1]),
                     .packages = "popbio") %dopar% {
                       # change to orginal code - sim.out is now a list, every rep is
                       #   its own entry. removed .combine = cbind from above
                       
                       # start pva run
                       #init.N <- simul.in[s, 1]
                       # shift into function argument
                       default <- ROI.PVA(n.reps,
                                      n.years,
                                      init.N      = ROI.simul.in[s, 1],
                                      ROI.K, 
                                      catastrophe.p,
                                      catastrophe.i,
                                      F.loop      = ROI.simul.in[s, 4],
                                      S.jw.loop   = ROI.simul.in[s, 3],
                                      S.ad.loop   = ROI.simul.in[s, 2],
                                      philopatry  = ROI.simul.in[s, 5],
                                      breed.fail)
                       
                     }
  
)

# Stop parallel processing
stopCluster(cl)

roi.names <- rep(as.character(ROI.pop.1985), (104652/9))
names(ROI.sim.out) <- roi.names


saveRDS(ROI.sim.out,
        "data/ROI_base_mod_1985_present_sims.RDS")
# version 2 of ROI 1985 base saved with breed.fail 0.5, ROI.K = 25000



#-----  ROI Plotting  -----
ROI.sim.out <- readRDS("data/ROI_base_mod_1985_present_sims_2.RDS")

roi.names <- rep(as.character(ROI.pop.1985), (104652/9))
names(ROI.sim.out) <- roi.names

years <- seq(1987, 2017, 1)

ROI.sim.out.5000 <- ROI.sim.out[names(ROI.sim.out) == "5000"]
ROI.sim.out.5000 <- do.call(cbind, ROI.sim.out.5000)
ROI.sim.out.5000.mean <- rowMeans(ROI.sim.out.5000)
ROI.sim.out.5000.lcl  <- apply(ROI.sim.out.5000, 1, quantile, probs = 0.1)
ROI.sim.out.5000.ucl  <- apply(ROI.sim.out.5000, 1, quantile, probs = 0.9)


ROI.sim.out.5500 <- ROI.sim.out[names(ROI.sim.out) == "5500"]
ROI.sim.out.5500 <- do.call(cbind, ROI.sim.out.5500)
ROI.sim.out.5500.mean <- rowMeans(ROI.sim.out.5500)
ROI.sim.out.5500.lcl  <- apply(ROI.sim.out.5500, 1, quantile, probs = 0.1)
ROI.sim.out.5500.ucl  <- apply(ROI.sim.out.5500, 1, quantile, probs = 0.9)

ROI.sim.out.6000 <- ROI.sim.out[names(ROI.sim.out) == "6000"]
ROI.sim.out.6000 <- do.call(cbind, ROI.sim.out.6000)
ROI.sim.out.6000.mean <- rowMeans(ROI.sim.out.6000)
ROI.sim.out.6000.lcl  <- apply(ROI.sim.out.6000, 1, quantile, probs = 0.1)
ROI.sim.out.6000.ucl  <- apply(ROI.sim.out.6000, 1, quantile, probs = 0.9)

ROI.sim.out.6500 <- ROI.sim.out[names(ROI.sim.out) == "6500"]
ROI.sim.out.6500 <- do.call(cbind, ROI.sim.out.6500)
ROI.sim.out.6500.mean <- rowMeans(ROI.sim.out.6500)
ROI.sim.out.6500.lcl  <- apply(ROI.sim.out.6500, 1, quantile, probs = 0.1)
ROI.sim.out.6500.ucl  <- apply(ROI.sim.out.6500, 1, quantile, probs = 0.9)

ROI.sim.out.7000 <- ROI.sim.out[names(ROI.sim.out) == "7000"]
ROI.sim.out.7000 <- do.call(cbind, ROI.sim.out.7000)
ROI.sim.out.7000.mean <- rowMeans(ROI.sim.out.7000)
ROI.sim.out.7000.lcl  <- apply(ROI.sim.out.7000, 1, quantile, probs = 0.1)
ROI.sim.out.7000.ucl  <- apply(ROI.sim.out.7000, 1, quantile, probs = 0.9)

ROI.sim.out.7500 <- ROI.sim.out[names(ROI.sim.out) == "7500"]
ROI.sim.out.7500 <- do.call(cbind, ROI.sim.out.7500)
ROI.sim.out.7500.mean <- rowMeans(ROI.sim.out.7500)
ROI.sim.out.7500.lcl  <- apply(ROI.sim.out.7500, 1, quantile, probs = 0.1)
ROI.sim.out.7500.ucl  <- apply(ROI.sim.out.7500, 1, quantile, probs = 0.9)

ROI.sim.out.8000 <- ROI.sim.out[names(ROI.sim.out) == "8000"]
ROI.sim.out.8000 <- do.call(cbind, ROI.sim.out.8000)
ROI.sim.out.8000.mean <- rowMeans(ROI.sim.out.8000)
ROI.sim.out.8000.lcl  <- apply(ROI.sim.out.8000, 1, quantile, probs = 0.1)
ROI.sim.out.8000.ucl  <- apply(ROI.sim.out.8000, 1, quantile, probs = 0.9)

ROI.sim.out.8500 <- ROI.sim.out[names(ROI.sim.out) == "8500"]
ROI.sim.out.8500 <- do.call(cbind, ROI.sim.out.8500)
ROI.sim.out.8500.mean <- rowMeans(ROI.sim.out.8500)
ROI.sim.out.8500.lcl  <- apply(ROI.sim.out.8500, 1, quantile, probs = 0.1)
ROI.sim.out.8500.ucl  <- apply(ROI.sim.out.8500, 1, quantile, probs = 0.9)

ROI.sim.out.9000 <- ROI.sim.out[names(ROI.sim.out) == "9000"]
ROI.sim.out.9000 <- do.call(cbind, ROI.sim.out.9000)
ROI.sim.out.9000.mean <- rowMeans(ROI.sim.out.9000)
ROI.sim.out.9000.lcl  <- apply(ROI.sim.out.9000, 1, quantile, probs = 0.1)
ROI.sim.out.9000.ucl  <- apply(ROI.sim.out.9000, 1, quantile, probs = 0.9)



# Combine into tibble


ROI.sim.out.df.long <- data.frame(years = rep(years, times = 9),
                             initN = rep(ROI.pop.1985, each = 31),
                             mean = c(ROI.sim.out.5000.mean,
                                      ROI.sim.out.5500.mean,
                                      ROI.sim.out.6000.mean,
                                      ROI.sim.out.6500.mean,
                                      ROI.sim.out.7000.mean,
                                      ROI.sim.out.7500.mean,
                                      ROI.sim.out.8000.mean,
                                      ROI.sim.out.8500.mean,
                                      ROI.sim.out.9000.mean),
                             lcl = c(ROI.sim.out.5000.lcl,
                                     ROI.sim.out.5500.lcl,
                                     ROI.sim.out.6000.lcl,
                                     ROI.sim.out.6500.lcl,
                                     ROI.sim.out.7000.lcl,
                                     ROI.sim.out.7500.lcl,
                                     ROI.sim.out.8000.lcl,
                                     ROI.sim.out.8500.lcl,
                                     ROI.sim.out.9000.lcl),
                             ucl = c(ROI.sim.out.5000.ucl,
                                     ROI.sim.out.5500.ucl,
                                     ROI.sim.out.6000.ucl,
                                     ROI.sim.out.6500.ucl,
                                     ROI.sim.out.7000.ucl,
                                     ROI.sim.out.7500.ucl,
                                     ROI.sim.out.8000.ucl,
                                     ROI.sim.out.8500.ucl,
                                     ROI.sim.out.9000.ucl)
)


write.table(ROI.sim.out.df.long,
            "data/ROI.sim.out.base.1985.means.csv",
            sep = ",",
            row.names = FALSE)

ROI.pop.ests <- data.frame(years = c(1987, 2017),
                          popsize = c((2*7000), (2*150)))


ROI.sim.out.plot <- ggplot(ROI.sim.out.df.long) +
  facet_wrap(~initN) + 
  geom_line(aes(x = years,
                y = mean,
                colour = factor(initN)))

ROI.sim.out.plot.ribbon <- ROI.sim.out.plot + 
  geom_ribbon(aes(x = years,
                  ymin = lcl,
                  ymax = ucl),
              alpha = 0.3) + 
  geom_point(data = ROI.pop.ests,
             aes(x = years,
                 y = popsize)) +
  xlab("Years") + 
  ylab("Curlew Population Size (individuals)") +
  ylim(0, 22000) + 
  theme(legend.position = "none",
        panel.background = element_rect(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())




#----- Obtaining predictions that predict less than 500 curlew  -----
bettermatch.5000 <- which(ROI.sim.out.5000[31 ,] < 200,
                     arr.ind = TRUE)
ROI.sim.out.5000.match <- ROI.sim.out.5000[, bettermatch.5000]
ROI.sim.out.5000.match.mean <- rowMeans(ROI.sim.out.5000.match)
ROI.sim.out.5000.match.lcl <- apply(ROI.sim.out.5000.match, 1,
                                   quantile, probs = 0.1)
ROI.sim.out.5000.match.ucl <- apply(ROI.sim.out.5000.match, 1,
                                   quantile, probs = 0.9)

bettermatch.5500 <- which(ROI.sim.out.5500[31 ,] < 200,
                          arr.ind = TRUE)
ROI.sim.out.5500.match <- ROI.sim.out.5500[, bettermatch.5500]
ROI.sim.out.5500.match.mean <- rowMeans(ROI.sim.out.5500.match)
ROI.sim.out.5500.match.lcl <- apply(ROI.sim.out.5500.match, 1,
                                    quantile, probs = 0.1)
ROI.sim.out.5500.match.ucl <- apply(ROI.sim.out.5500.match, 1,
                                    quantile, probs = 0.9)

bettermatch.6000 <- which(ROI.sim.out.6000[31 ,] < 200,
                          arr.ind = TRUE)
ROI.sim.out.6000.match <- ROI.sim.out.6000[, bettermatch.6000]
ROI.sim.out.6000.match.mean <- rowMeans(ROI.sim.out.6000.match)
ROI.sim.out.6000.match.lcl <- apply(ROI.sim.out.6000.match, 1,
                                    quantile, probs = 0.1)
ROI.sim.out.6000.match.ucl <- apply(ROI.sim.out.6000.match, 1,
                                    quantile, probs = 0.9)

bettermatch.6500 <- which(ROI.sim.out.6500[31 ,] < 200,
                          arr.ind = TRUE)
ROI.sim.out.6500.match <- ROI.sim.out.6500[, bettermatch.6500]
ROI.sim.out.6500.match.mean <- rowMeans(ROI.sim.out.6500.match)
ROI.sim.out.6500.match.lcl <- apply(ROI.sim.out.6500.match, 1,
                                    quantile, probs = 0.1)
ROI.sim.out.6500.match.ucl <- apply(ROI.sim.out.6500.match, 1,
                                    quantile, probs = 0.9)

bettermatch.7000 <- which(ROI.sim.out.7000[31 ,] < 200,
                          arr.ind = TRUE)
ROI.sim.out.7000.match <- ROI.sim.out.7000[, bettermatch.7000]
ROI.sim.out.7000.match.mean <- rowMeans(ROI.sim.out.7000.match)
ROI.sim.out.7000.match.lcl <- apply(ROI.sim.out.7000.match, 1,
                                    quantile, probs = 0.1)
ROI.sim.out.7000.match.ucl <- apply(ROI.sim.out.7000.match, 1,
                                    quantile, probs = 0.9)

bettermatch.7500 <- which(ROI.sim.out.7500[31 ,] < 200,
                          arr.ind = TRUE)
ROI.sim.out.7500.match <- ROI.sim.out.7500[, bettermatch.7500]
ROI.sim.out.7500.match.mean <- rowMeans(ROI.sim.out.7500.match)
ROI.sim.out.7500.match.lcl <- apply(ROI.sim.out.7500.match, 1,
                                    quantile, probs = 0.1)
ROI.sim.out.7500.match.ucl <- apply(ROI.sim.out.7500.match, 1,
                                    quantile, probs = 0.9)

bettermatch.8000 <- which(ROI.sim.out.8000[31 ,] < 200,
                          arr.ind = TRUE)
ROI.sim.out.8000.match <- ROI.sim.out.8000[, bettermatch.8000]
ROI.sim.out.8000.match.mean <- rowMeans(ROI.sim.out.8000.match)
ROI.sim.out.8000.match.lcl <- apply(ROI.sim.out.8000.match, 1,
                                    quantile, probs = 0.1)
ROI.sim.out.8000.match.ucl <- apply(ROI.sim.out.8000.match, 1,
                                    quantile, probs = 0.9)

bettermatch.8500 <- which(ROI.sim.out.8500[31 ,] < 200,
                          arr.ind = TRUE)
ROI.sim.out.8500.match <- ROI.sim.out.8500[, bettermatch.8500]
ROI.sim.out.8500.match.mean <- rowMeans(ROI.sim.out.8500.match)
ROI.sim.out.8500.match.lcl <- apply(ROI.sim.out.8500.match, 1,
                                    quantile, probs = 0.1)
ROI.sim.out.8500.match.ucl <- apply(ROI.sim.out.8500.match, 1,
                                    quantile, probs = 0.9)

bettermatch.9000 <- which(ROI.sim.out.9000[31 ,] < 200,
                          arr.ind = TRUE)
ROI.sim.out.9000.match <- ROI.sim.out.9000[, bettermatch.9000]
ROI.sim.out.9000.match.mean <- rowMeans(ROI.sim.out.9000.match)
ROI.sim.out.9000.match.lcl <- apply(ROI.sim.out.9000.match, 1,
                                    quantile, probs = 0.1)
ROI.sim.out.9000.match.ucl <- apply(ROI.sim.out.9000.match, 1,
                                    quantile, probs = 0.9)



ROI.sim.out.df.match <- data.frame(years = rep(years, times = 9),
                                  initN = rep(ROI.pop.1985, each = 31),
                                  mean = c(ROI.sim.out.5000.match.mean,
                                           ROI.sim.out.5500.match.mean,
                                           ROI.sim.out.6000.match.mean,
                                           ROI.sim.out.6500.match.mean,
                                           ROI.sim.out.7000.match.mean,
                                           ROI.sim.out.7500.match.mean,
                                           ROI.sim.out.8000.match.mean,
                                           ROI.sim.out.8500.match.mean,
                                           ROI.sim.out.9000.match.mean),
                                  lcl = c(ROI.sim.out.5000.match.lcl,
                                          ROI.sim.out.5500.match.lcl,
                                          ROI.sim.out.6000.match.lcl,
                                          ROI.sim.out.6500.match.lcl,
                                          ROI.sim.out.7000.match.lcl,
                                          ROI.sim.out.7500.match.lcl,
                                          ROI.sim.out.8000.match.lcl,
                                          ROI.sim.out.8500.match.lcl,
                                          ROI.sim.out.9000.match.lcl),
                                  ucl = c(ROI.sim.out.5000.match.ucl,
                                          ROI.sim.out.5500.match.ucl,
                                          ROI.sim.out.6000.match.ucl,
                                          ROI.sim.out.6500.match.ucl,
                                          ROI.sim.out.7000.match.ucl,
                                          ROI.sim.out.7500.match.ucl,
                                          ROI.sim.out.8000.match.ucl,
                                          ROI.sim.out.8500.match.ucl,
                                          ROI.sim.out.9000.match.ucl)
)

write.table(ROI.sim.out.df.match,
            "data/ROI_sims_1985_2017_lowest_v1.csv",
            sep = ",",
            row.names = FALSE)

ROI.sim.out.match.plot <- ggplot(ROI.sim.out.df.match) +
  facet_wrap(~initN) + 
  geom_line(aes(x = years,
                y = mean,
                colour = factor(initN)))

ROI.sim.out.match.plot.ribbon <- ROI.sim.out.match.plot + 
  geom_ribbon(aes(x = years,
                  ymin = lcl,
                  ymax = ucl),
              alpha = 0.3) + 
  geom_point(data = ROI.pop.ests,
             aes(x = years,
                 y = popsize)) +
  xlab("Years") + 
  ylab("Curlew Population Size (individuals)") +
  ylim(0, 20000) + 
  theme(legend.position = "none",
        panel.background = element_rect(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())



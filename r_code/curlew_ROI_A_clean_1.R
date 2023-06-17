##### Head-starting Curlew Population Viability Analysis  #####

# Head-Starting Scenario A - courtesy of Becs Lee
# Underlying productivity does not change, low productivity continues
#   in the long-term
#   1. Scenario A1: no head-starting - this is the BASELINE MODEL
#       - 0 clutches added
#   2. Scenario A2: enough head-starting to prevent decline (i.e. minimum 
#       release numbers), duration 5 years, HS birds survive & return
#       in line with wild-reared curlew
#       - this was tricky, as we were just visualizing from the figures
#       - looking for a flattened line for the first 8 years
#         as this was under the assumption curlew don't start breeding
#         regularly until they are 3 years in age - which may be a wrong 
#         assumption
#   3. Scenario A3: maximum head-starting that is reasonably possible
#       (e.g. collecting 20 clutches), duration 5 years, HS birds
#       survive & return in line with wild-reared curlew
#       For NI this was 60
#   4. Scenario A4: worst case scenario - maximum collection strategy is taken
#       (same as A3, 20 clutches), duration 5 years, HS birds do not survive
#       or do not return (make no contribution to population)
#

# The basis of the main model adapted from Shoemaker and Oppel.
# Model details:
# Demographic vital rates were extracted from the literature and 
#   agreed upon in discussion with an RSPB and colleagues steering group.
# These vital rates will be used to create a table of all possible combinations
#   each of which will be simulated over to predict 
# Initial population, N, will be the current estimates (2021) of curlew across
#   NI & ROI. Depending on spatial scale, these values can be restricted to 
#   region, local site, all island.
# 



#-----  Load Libraries  -----
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


#-----  Parameters  -----

# Population Size:
#   - the population size is varied for each sitewe are modelling
#   - The initial NI pop size is esimtated to be between 100 - 250.
#   - the variable name changes when modelling all of NI (NI.pop.size)
#     and at site level (site.pop)
#   - populations for sites:
#         - glenwherry = 57
#         - lle        = 37
#

ROI.pop.size <- seq(100, 250, 50)


# Survival Rates & Curlew Related Parameters
# All sources for paramaters can be found in the appendix for the report
# or the separate word document

S.ad            <- seq(0.78, 0.92, 0.04)    # survival adult

S.jw            <- seq(0.20, 0.40, 0.05)    # survival of wild fledgling

S.jhs           <- seq(0.20, 0.40, 0.05)    # survival of head-started fledge

S.sub           <- seq(0.78, 0.92, 0.04)    # survival subadult

S.imm           <- seq(0.40, 0.70, 0.10)    # survival of immature

prod.rates      <- seq(0.14, 0.52, 0.06)    # productivity rate

hs.clutches     <- 0                        # num of clutches to add
hs.pop          <- (hs.clutches*4)          # num of birds added -
# based on est. 4 eggs per clutch
# 0 eggs is base model i.e. current
# state
# Max clutches (scen 3) = 60

capt.fail       <- 0.2                      # % birds collected that are
# are not added to the pop due to
# captive breeding failure etc

intro.years     <- 5                        # yrs in which HS birds are 
# introduced

phil.p           <- 0.4                     # set single value to speed up model
#phil.p          <- seq(0.2, 0.4, 0.6)       # philopatry - probability of 
# return for first breeding
# applies to first & second years

return.p        <- 0.9                      # probability breeding birds return
#return.p       <- seq(0.7, 0.9, 0.1)

# Environmental Parameters
K               <- 10000                    # carrying capacity

site.K          <- 10000                    # carrying capacity of site

# K <- site.K                               # replace K with site.K so you
# don't have to change model code

# Stochastic Parameters
SD.lambda       <- 0.1                      # SD of growth rate
# enables simulation of 
# stochasticity in lambda by
# randomly drawing a value
# from a normal dist between 0 &
# lambda.max, with variation of SD

catastrophe.p   <- 0.001                    # prob of catastrophic event

catastrophe.i   <- 0.25                     # impact of catastrophe on pop
# essentially reduces by 1/4
# multiply pop by 0.25

breed.fail      <- 0.4                      # prob breeding season is total 
# failure

extinct.cutoff  <- c(1, 10, 30, 60)

acceptable.risk <- 0.1



# PVA settings
n.years         <- 30                       # number of years to predict pop
n.reps          <- 500                      # number of simulations for each 
# combo of params


#-----  Defining Functions for PVA  -----

# Stage Matrix
#   the structure of the population model with the different stages
#   fledgling = survival of fledgling
#   jw        = wild chicks age fledge - 1
#   jhs       = head-started chicks age fledge - 1
#   immature  = age 1 - 2
#   subadult  = age 2 - 3
#   adult     = age 3 +
# Several reports suggest curlew can breed from age 1 + years
# To account for this, we keep the productivity the same across all age groups
#   but multiply that productivity by a proportion of the pop within that
#   age group - so not all birds in that group breed.
#     - imm set to 0.5% of age group breed
#     - sub
# 


bird.matrix <- expression(
  # f_w       f_hs      imm                     sub                     adult  
  0,          0,    (prod.rates*0.5)*0.005, (prod.rates*0.5)*0.1, (prod.rates*0.5),
  0,          0,     0,                       0,                      0,
  S.jw*phil.p, (S.jhs*(1-capt.fail)),  0,  0,                      0,
  0,          0,      S.imm*return.p,         0,                      0,
  0,          0,      0,                      S.sub*return.p,     S.ad*return.p
)


# This matrix is for scenario 4 in the model - where
#   there is reduced relay probability and no returning birds 
# Bird matrix with reduced productivity for relay
# Relay potential reduced by 50%
bird.matrix.2 <- expression(
  # f_w       f_hs      imm                     sub                     adult  
  0,          0,    (prod.rates*0.5)*0.005, (prod.rates*0.5)*0.1, (prod.rates*0.5)*(1 - (((init.N - hs.clutches)*0.5)/100)),
  0,          0,     0,                       0,                      0,
  (S.jw*phil.p)*0, (S.jhs*(1-capt.fail)),  0,  0,                      0,
  0,          0,      S.imm*return.p,         0,                      0,
  0,          0,      0,                      S.sub*return.p,     S.ad*return.p
)



# Density dependent Ricker Model to simulate population growth

Ricker <- function(prev.abund) {
  
  # function for calculating next year abundance
  # includes environmental/demographic stochasticity
  
  prev.abund * exp(log(rnorm(1,
                             R.max,
                             SD.lambda))
                   * (1 - (prev.abund/K)))
  
}


R.max <- 1



#----- Create Combinations of Demographic Parameters  -----
#
# Use expand grid to create a table of all possible parameter combinations
# Initially, this resulted in a ridiculous number of combinations due to
# the range in vales for each parameter.
# Speaking to Steffen and Gillian, we whittled the ranges down for some,
#   but this can be further reduced in the param set-up above. Just limit the
#   number of values for each paramater

simul.in  <- expand.grid(init.N       = ROI.pop.size,
                         S.ad.loop    = S.ad,
                         S.jw.loop    = S.jw,
                         S.jhs.loop   = S.jhs,
                         S.imm.loop   = S.imm,
                         S.sub.loop   = S.sub,
                         F.loop       = prod.rates,
                         philopatry   = phil.p,
                         sitefidelity = return.p,
                         hs.numbers   = hs.pop,
                         hs.years     = intro.years)

names(simul.in) <- c("initial.pop",
                     "S.ad",
                     "S.jw",
                     "S.jhs",
                     "S.imm",
                     "S.sub",
                     "Prod.Rates",
                     "philopatry",
                     "site.fidelity",
                     "number.HS.birds",
                     "years.HS")

# Save to file - change names to whatever you need
write.table(simul.in,
            "data/curlew_ROI_A_params_0c.csv",
            sep = ",",
            row.names = FALSE)



# Simplifying the parameter combinations:
#   After discussions with the strategy group, it was decided to run a more
#   simplified model and reduce the number of parameter combinations further
#   by only looking at the maximum, minimum and mean values for each.
#   This is how I extracted them from the above table
simul.in.min <- simul.in %>%
  summarise_if(is.numeric, min) 
simul.in.mean <- simul.in %>%
  summarise_if(is.numeric, mean) 
simul.in.max <- simul.in %>%
  summarise_if(is.numeric, max)

# Combine simplified params into a data frame
simul.in.simple <- as.data.frame(
  rbind(simul.in.min,
        simul.in.mean,
        simul.in.max)
)

names(simul.in.simple) <- c("initial.pop",
                            "S.ad",
                            "S.jw",
                            "S.jhs",
                            "S.imm",
                            "S.sub",
                            "Prod.Rates",
                            "philopatry",
                            "site.fidelity",
                            "number.HS.birds",
                            "years.HS")

# Save to file - change names to whatever you need
write.table(simul.in.simple,
            "data/curlew_ROI_A_simple_params_0c.csv",
            sep = ",",
            row.names = FALSE)



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
  S.ad.loop,
  S.jw.loop,
  S.jhs.loop,
  S.imm.loop,
  S.sub.loop,
  philopatry,
  sitefidelity,
  hs.numbers,
  hs.years,
  capt.fail,
  breed.fail
) {
  
  # create array to store pop values
  # rows = yars, cols = reps
  Pop.Array <- array(0,
                     dim = c((n.years + 1),
                             n.reps))
  
  # start loop of reps
  for (rep in 1:n.reps) {
    
    # initial abundance + HS birds
    # don't add hs birds as pop count is only breeding birds
    Pop.Array[1, rep]   <- (init.N*2) 
    
    
    # Loop through the years
    for (y in 2:(n.years + 1)) {
      
      # create matrix with vital rates
      # Check if ruined breeding season (20%) - varies every sim
      F.Year <- ifelse(rbinom(1, 1, breed.fail) == 1,
                       F.loop * 0.2,
                       F.loop)
      
      # Create list of vital rates
      bird.vr <- list(init.N = init.N,
                      prod.rates  = F.Year,
                      S.ad        = S.ad.loop,
                      S.jw        = S.jw.loop,
                      S.jhs       = S.jhs.loop,
                      S.imm       = S.imm.loop,
                      S.sub       = S.sub.loop,
                      phil.p      = philopatry,
                      return.p    = sitefidelity,
                      capt.fail   = capt.fail,
                      hs.clutches = hs.clutches
      )
      
      # Create leslie matrix
      if (y <= hs.years) {
        a <- matrix(sapply(bird.matrix, # change to suit model scenarios
                           # for example, use bird.matrix.2 for 
                           # scenario A4
                           eval,
                           bird.vr,
                           NULL),
                    nrow = sqrt(length(bird.matrix)),
                    byrow = TRUE)
      } else {
        a <- matrix(sapply(bird.matrix,
                           eval,
                           bird.vr,
                           NULL),
                    nrow = sqrt(length(bird.matrix)),
                    byrow = TRUE)
      }
      
      # Starting pop vector
      # 0 wild chicks, as pop counts only report adults
      # introduce HS chicks
      # After first year, use calculated vectors from pop.projections
      # If only counting breeding pairs, only looking for adults
      # for scenario 4 multiply all f_w & f_hs by 0 every year
      if (y == 2) {
        pop.size  <- c((F.loop * 0.5 * (init.N * 2)),                       #f_w
                       (hs.numbers * (1-capt.fail)),                        #f_hs
                       ((hs.numbers * (1-capt.fail)) * S.jhs.loop * philopatry + 
                          (F.loop * 0.5 * (init.N *2) * philopatry)), #imm
                       0,                                              # sub
                       ((init.N * 2) * S.ad.loop * sitefidelity))      #adu
      } else if (y == 3) {
        pop.size  <- c(projections$stage.vectors[1, y-1],
                       (hs.numbers * (1 - capt.fail)),
                       projections$stage.vectors[3, y-1],
                       projections$stage.vectors[4, y-1],
                       projections$stage.vectors[5, y-1])
      } else if (y == 4) {
        pop.size  <- c(projections$stage.vectors[1, y-1],
                       (hs.numbers * (1 - capt.fail)),
                       projections$stage.vectors[3, y-1],
                       projections$stage.vectors[4, y-1],
                       projections$stage.vectors[5, y-1])
      } else if (y == 5) {
        pop.size  <- c(projections$stage.vectors[1, y-1],
                       (hs.numbers * (1 - capt.fail)),
                       projections$stage.vectors[3, y-1],
                       projections$stage.vectors[4, y-1],
                       projections$stage.vectors[5, y-1])
      } else {
        pop.size  <- c(projections$stage.vectors[, y-1])
      }
      
      
      
      # Get projections for years
      projections <- pop.projection(a,
                                    n           = pop.size,
                                    iterations  = 30)
      
      # Stochastic Ricker Model
      R.max <- projections$lambda
      
      # Calculate abundance
      next.year   <- max(0,
                         trunc(Ricker(Pop.Array[y - 1, rep])))
      
      # Catastrophe
      if(runif(1) < catastrophe.p) next.year <- next.year * catastrophe.i
      
      # Put in array
      Pop.Array[y, rep] <- next.year
      
    } # end year
    
  }   # end rep
  
  return(Pop.Array)
  
}



#-----  Run Sims with Parallel Processing -----

# Boost speed using parallel processing
cl <- makeCluster(8)
registerDoParallel(cl, cores = 8)

# Start loop
system.time(
  
  ROI.sim.out <- foreach(s = c(1:dim(simul.in.simple)[1]),
                        .packages = "popbio") %dopar% {
                          
                          # start PVA run
                          init.N           = simul.in.simple[s, 1]
                          default <- PVA(
                            n.reps,
                            n.years,
                            init.N,
                            K,
                            catastrophe.p,
                            catastrophe.i,
                            F.loop           = simul.in.simple[s, 7],
                            S.ad.loop        = simul.in.simple[s, 2],
                            S.jw.loop        = simul.in.simple[s, 3],
                            S.jhs.loop       = simul.in.simple[s, 4],
                            S.imm.loop       = simul.in.simple[s, 5],
                            S.sub.loop       = simul.in.simple[s, 6],
                            philopatry       = simul.in.simple[s, 8],
                            sitefidelity     = simul.in.simple[s, 9],
                            hs.numbers       = simul.in.simple[s, 10],
                            hs.years         = simul.in.simple[s, 11],
                            capt.fail,
                            breed.fail
                          )
                          
                        }
  
)


# Stop processing
stopCluster(cl)



###########

#-----  Reduced number of runs to a single run  -----
# AKA the simplified model
n.reps <- 1
s = 1, 2, 3 # select a row from the table of parameters for this single run
# when running the simplified model, there are only 3 rows

testPVA <- PVA(
  n.reps,
  n.years,
  init.N = simul.in.simple[s, 1],             # when running full model 
  # change to simul.in[s, 1]
  K,
  catastrophe.p,
  catastrophe.i,
  F.loop           = simul.in.simple[s, 7],
  S.ad.loop        = simul.in.simple[s, 2],
  S.jw.loop        = simul.in.simple[s, 3],
  S.jhs.loop       = simul.in.simple[s, 4],
  S.imm.loop       = simul.in.simple[s, 5],
  S.sub.loop       = simul.in.simple[s, 6],
  philopatry       = simul.in.simple[s, 8],
  sitefidelity     = simul.in.simple[s, 9],
  hs.numbers       = simul.in.simple[s, 10],
  hs.years         = simul.in.simple[s, 11],
  capt.fail,
  breed.fail
)


# When using the simplified model paramaters - we classify these as
#   min, max, mean
# Write out testPVA to one of the above

min.PVA <- testPVA   # where s = 1
med.PVA <- testPVA   # where s = 2
max.PVA <- testPVA   # where s = 3

# create vector of years
years <- seq(2023, 2053, 1)

PVA_ROI_base_model <- as.data.frame(
  cbind(
    years,
    min.PVA,
    med.PVA,
    max.PVA
  )
)


# Write to file
# Change scenario name in file and model
write.table(PVA_ROI_base_model,
            "data/chs_ROI_scenario_A_PVA_base.csv",
            sep = ",",
            row.names = FALSE)



# I couldn't figure a function to combine all of the post model run bits. It
#   got even worse for some of the ROI stuff when I had to flip the tables
#   into long form for plotting.
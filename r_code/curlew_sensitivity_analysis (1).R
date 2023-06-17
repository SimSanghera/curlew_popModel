#----- Sensitivity & Elasticity Analysis  -----
# Use sensitivyt and elasticity function sin popbio package to look at the 
#  sensitivity of each factor in the age matrix
#  - survival, prods, phil, return
# Create age matrices with different values for each param
# Create sensitivity matrices for each & elasticity
# Compare differences between and plot




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

# Population Size
# Population Size
#NI.pop.size   <- seq(100, 250, 50)      # breeding pairs
site.pop      <- 150                     # number of breed pairs on site
# glen  = 57
# lle = 37
# low = 4
# lowermed = 8
# med = 12
# uppermed  = 16
# high = 20

#ROI.pop.size  <- seq(100, 120, 10)


# Survival Rates & Curlew Related Parameters
S.ad          <- seq(0.78, 0.92, 0.04)  # Cook et al., 2021
# Roodbergen et al., 2021

S.jw          <- seq(0.20, 0.40, 0.05)  # Wild bird survival
# Cook et al., 2021

S.jhs         <- seq(0.20, 0.40, 0.05)  # HS bird survival
# Assume same as wild for now
# Geoff Hilton (pers comms)

S.sub         <- seq(0.78, 0.92, 0.04)  # subadult survival

S.imm         <- seq(0.40, 0.70, 0.10)  # survival of immature (age fledge - 1)

prod.rates    <- seq(0.14, 0.52, 0.06)  # Grant et al., 1999
# Harris et al., 2019
# Zielonka et al., 2019
# Irish Breeding Bird Report

hs.clutches   <- 20
hs.pop        <- (4*hs.clutches)        # baseline - 0 HS
# capt fail = 20% - start with 10 clutches then increase

capt.fail     <- 0.2                    # The number of HS birds that
# fail during rearing & can't be
# released

intro.years   <- 5           # years in which hs birds will be
# introduced

phil.p        <- 0.5                    # probability of first years return
# to release site
# philopatry is still possible in wild
# bred birds, but not HS

return.p      <- 0.95                   # probability of breeding birds 
# returning to release site


# Environmental Parameters
site.K          <- 1000                  # carrying capacity of Glenwherry
# other sites will have lower K
low.k           <- 5 #10, 15, 20, 50, 100

K <- site.K

# Stochastic Parameters
SD.lambda     <- 0.10                   # Using SD of growth rate enables 
# simulation of stochastic lambda
# by randomly drawing a value from 
# a normal distribution between 0
# and lambda max, with a variation
# of the SD.

catastrophe.p <- 0.001                  # probability of a catastrophic event
# occurring - flood, earthquake, 
# human-disturbance

catastrophe.i <- 0.25                   # the impact of the catastrophe on
# the pop. multiply pop by 0.25,
# essentially leaving 1/4 of total
# pop alive

breed.fail    <- 0.2                    # prob the breeding season is a 
# total failure within the pop.
# This could potentially be 
# calculated if annual monitoring 
# records are available - eg. Glenwherry

extinct.cutoff  <- c(1, 10, 30, 60)      # below this number, pop is extinct

acceptable.risk <- 0.1                  # risk considered acceptable to 
# managers

# sex.ratio     <- not viable for this model. assume 50/50 sex ratio


# PVA Settings
n.years         <- 30                   # number of years to simulate
n.reps          <- 500                  # number of simulations for each 
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
  S.jw*phil.p, (S.jhs*capt.fail)*phil.p,  0,  0,                      0,
  0,          0,      S.imm*return.p,         0,                      0,
  0,          0,      0,                      S.sub*return.p,     S.ad*return.p
)


# Bird matrix with reduced productivity for relay
# Relay potential reduced by 50%
# Remove all wild births - f_w
# Remove *0.5 for prod.rates as this is for single bird rates
bird.matrix.2 <- expression(
  # f_w       f_hs      imm                     sub                     adult  
  0,          0,    (prod.rates*0.5)*0.005, (prod.rates*0.5)*0.1, (prod.rates*0.5)*(1 - (((init.N - hs.clutches)*0.5)/100)),
  0,          0,     0,                       0,                      0,
  (S.jw*phil.p)*0, (S.jhs*capt.fail)*phil.p,  0,  0,                      0,
  0,          0,      S.imm*return.p,         0,                      0,
  0,          0,      0,                      S.sub*return.p,     S.ad*return.p
)


# To calculate correct relaying proportion:
# - subtract number of clutches from breeding pop 
#   e.g. 57 breeding pairs, take away 25 clutches
#       32 pairs single lay = 100%
#       25 pairs 50% relay = 12.5 birds relay
#       100% - (x*50%)
#       x = site.pop-no.clutches
# (1 - ((site.pop - hs.clutches)*0.5) - eqn for 

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


# Calculate proportion of simulations wehre species go extinct
extinction.bysim <- function(simdata,
                             threshold) {
  
  # threshold number set above
  sum(apply(default,
            2,
            function(t) min(t) < threshold)) / ncol(simdata)
  
}


#-----  Create Combinations of Demographic Parameters -----

simul.in  <- expand.grid(init.N       = site.pop,
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

simul.in.min <- simul.in %>%
  summarise_if(is.numeric, min) 
simul.in.mean <- simul.in %>%
  summarise_if(is.numeric, mean) 
simul.in.max <- simul.in %>%
  summarise_if(is.numeric, max)

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

write.table(simul.in.simple,
            "data/curlew_HS_ROI_simple_params_A1.csv",
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
        a <- matrix(sapply(bird.matrix.2,   # reduced relay & no return
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

rep = 1
s = 1
y = 2
init.N           = simul.in.simple[s, 1]
F.loop           = simul.in.simple[s, 7]
S.ad.loop        = simul.in.simple[s, 2]
S.jw.loop        = simul.in.simple[s, 3]
S.jhs.loop       = simul.in.simple[s, 4]
S.imm.loop       = simul.in.simple[s, 5]
S.sub.loop       = simul.in.simple[s, 6]
philopatry       = simul.in.simple[s, 8]
sitefidelity     = simul.in.simple[s, 9]
hs.numbers       = simul.in.simple[s, 10]
hs.years         = simul.in.simple[s, 11]


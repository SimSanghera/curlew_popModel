######  Head-Starting Curlew Population Viability Analysis #####

# Head-Starting Scenario A - courtesy of Becs Lee
# Underlying productivity does not change, low productivity continues
#   in the long-term
#   1. Scenario A1: no head-starting
#       Baseline model. Upon the advice of SO, changed model params to reduce
#       number of values for a) computational efficiency and b) irrelevance
#       to model output.
#       The baseline model will be split into 6 different models:
#         a - low philopatry, high return
#         b - medium phil, high return
#         c - high phil, high return
#         d - low phil, low return
#         e - med phil, low return
#         f - high phil, low return
#         Low return is still a high prob at 0.7

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




#-----  Parameters  -----

# Population Size
# Population Size
NI.pop.size   <- seq(150, 250, 50)      # breeding pairs

ROI.pop.size  <- 110


# Survival Rates & Curlew Related Parameters
S.ad          <- seq(0.78, 0.92, 0.04)  # Cook et al., 2021
# Roodbergen et al., 2021

S.jw          <- seq(0.20, 0.40, 0.05)  # Wild bird survival
# Cook et al., 2021

S.jhs         <- seq(0.20, 0.40, 0.05)  # HS bird survival
# Assume same as wild for now
# Geoff Hilton (pers comms)

S.sub         <- seq(0.78, 0.92, 0.04)  # subadult survival

S.imm         <- seq(0.40, 0.70, 0.15)  # survival of immature (age fledge - 1)

prod.rates    <- seq(0.14, 0.52, 0.06)  # Grant et al., 1999
# Harris et al., 2019
# Zielonka et al., 2019
# Irish Breeding Bird Report


hs.clutches   <- 40
hs.pop        <- (4* hs.clutches)                   # baseline - 0 HS

capt.fail     <- 0.2                    # The number of HS birds that
# fail during rearing & can't be
# released

intro.years   <- 5           # years in which hs birds will be
# introduced

phil.p        <- 0.5                    # probability of first years return
# to release site

return.p      <- 0.90                   # probability of breeding birds 
# returning to release site


# Environmental Parameters
NI.K          <- 10000                  # carrying capacity of NI

ROI.K         <- 10000                  # carrying capacity of ROI

K <- NI.K

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

extinct.cutoff  <- c(1, 10, 30, 60)     # below this number, pop is extinct
# See notes in Joplin

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
  0,          0,    (prod.rates*0.5)*0.005, (prod.rates*0.5)*0.1, prod.rates*0.5,
  0,          0,     0,                       0,                      0,
  S.jw*phil.p, (S.jhs*phil.p),  0,  0,                      0,
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
                   * (1 - (prev.abund/NI.K)))
  
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
# hs years not needed in simul.in

simul.in  <- expand.grid(init.N       = NI.pop.size,
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
                     "hs.years")

write.table(simul.in,
            "data/NI_HS_A2_params_3.csv",
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
    Pop.Array[1, rep]   <- (init.N*2)
    
    
    # Loop through the years
    for (y in 2:(n.years + 1)) {
      
      # create matrix with vital rates
      # Check if ruined breeding season - varies every sim
      F.Year <- ifelse(rbinom(1, 1, breed.fail) == 1,
                       F.loop * 0.5,
                       F.loop)
      
      # Create list of vital rates
      bird.vr <- list(init.N      = init.N,
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
      a <- matrix(sapply(bird.matrix,
                         eval,
                         bird.vr,
                         NULL),
                  nrow = sqrt(length(bird.matrix)),
                  byrow = TRUE)
      
      # Starting pop vector
      # 0 wild chicks, as pop counts only report adults
      # introduce HS chicks
      # After first year, use calculated vectors from pop.projections
      if (y == 2) {
        pop.size  <- c((F.loop*0.5*(init.N*2)),                                # f_w
                       (hs.numbers * (1 - capt.fail)),                         # f_hs
                       ((hs.numbers*(1 - capt.fail)) * S.jhs.loop * philopatry), # imm
                       0,                                                      # sub
                       ((init.N*2)*S.ad.loop*sitefidelity))                    # ad
      } else if (y == 3) {
        pop.size  <- c(projections$stage.vectors[1, y-1],
                       (hs.numbers * (1 - capt.fail)),
                       projections$stage.vectors[3, y-1],
                       projections$stage.vectors[4, y-1],
                       projections$stage.vectors[5, y-1])
      } else if (y == 4) {
        pop.size  <- c(projections$stage.vectors[1, y-1],
                       (hs.numbers * (1- capt.fail)),
                       projections$stage.vectors[3, y-1],
                       projections$stage.vectors[4, y-1],
                       projections$stage.vectors[5, y-1])
      } else if (y == 5) {
        pop.size  <- c(projections$stage.vectors[1, y-1],
                       (hs.numbers * (1- capt.fail)),
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

# rep = 1
# y = 2
# s = 10
# init.N           = simul.in[1, 1]
# F.loop           = simul.in[s, 7]
# S.ad.loop        = simul.in[s, 2]
# S.jw.loop        = simul.in[s, 3]
# S.jhs.loop       = simul.in[s, 4]
# S.imm.loop       = simul.in[s, 5]
# S.sub.loop       = simul.in[s, 6]
# philopatry       = simul.in[s, 8]
# sitefidelity     = simul.in[s, 9]
# hs.numbers       = simul.in[s, 10]
# hs.years         = simul.in[s, 11]


#-----  Run Sims with Parallel Processing -----

# Boost speed using parallel processing
cl <- makeCluster(8)
registerDoParallel(cl, cores = 8)

# Start loop
system.time(
  
  NI.sim.out.35 <- foreach(s = c(1:dim(simul.in)[1]),
                     .packages = "popbio") %dopar% {
                       
                       # start PVA run
                       init.N           = simul.in[s, 1]
                       default <- PVA(
                         n.reps,
                         n.years,
                         init.N,
                         K,
                         catastrophe.p,
                         catastrophe.i,
                         F.loop           = simul.in[s, 7],
                         S.ad.loop        = simul.in[s, 2],
                         S.jw.loop        = simul.in[s, 3],
                         S.jhs.loop       = simul.in[s, 4],
                         S.imm.loop       = simul.in[s, 5],
                         S.sub.loop       = simul.in[s, 6],
                         philopatry       = simul.in[s, 8],
                         sitefidelity     = simul.in[s, 9],
                         hs.numbers       = simul.in[s, 10],
                         hs.years         = simul.in[s, 11],
                         capt.fail,
                         breed.fail
                       )
                       
                     }
  
)

# Stop processing
stopCluster(cl)

# What to name each element? intro bird numbers?


# Write to save file
saveRDS(NI.sim.out.35,
        "data/NI_A2_diffpops_35cl_1.RDS")



#----- Population plots -----
NI.sim.out.10 <- readRDS("data/NI_A2_diffpops_10cl_1.RDS")
NI.sim.out.15 <- readRDS("data/NI_A2_diffpops_15cl_1.RDS")


names.NI.sim.out <- rep(as.character(c(NI.pop.size)), (14175/3))

names(NI.sim.out.10) <- names.NI.sim.out
names(NI.sim.out.15) <- names.NI.sim.out
names(NI.sim.out.20) <- names.NI.sim.out
names(NI.sim.out.25) <- names.NI.sim.out
names(NI.sim.out.30) <- names.NI.sim.out
names(NI.sim.out.35) <- names.NI.sim.out


years <- seq(2023, 2053, 1)

# NI.sim.out.100 <- NI.sim.out[names(NI.sim.out) == "100"]
# NI.sim.out.100 <- do.call(cbind, NI.sim.out.100)
# NI.sim.out.100.mean <- rowMeans(NI.sim.out.100)
# NI.sim.out.100.lcl  <- apply(NI.sim.out.100, 1, quantile, probs = 0.1)
# NI.sim.out.100.ucl  <- apply(NI.sim.out.100, 1, quantile, probs = 0.9)

NI.sim.out.150.10 <- NI.sim.out[names(NI.sim.out.10) == "150"]
NI.sim.out.150.10 <- do.call(cbind, NI.sim.out.150.10)
NI.sim.out.150.10.mean <- rowMeans(NI.sim.out.150.10)
NI.sim.out.150.10.lcl  <- apply(NI.sim.out.150.10, 1, quantile, probs = 0.1)
NI.sim.out.150.10.ucl  <- apply(NI.sim.out.150.10, 1, quantile, probs = 0.9)

NI.sim.out.200.10 <- NI.sim.out[names(NI.sim.out.10) == "200"]
NI.sim.out.200.10 <- do.call(cbind, NI.sim.out.200.10)
NI.sim.out.200.10.mean <- rowMeans(NI.sim.out.200.10)
NI.sim.out.200.10.lcl  <- apply(NI.sim.out.200.10, 1, quantile, probs = 0.1)
NI.sim.out.200.10.ucl  <- apply(NI.sim.out.200.10, 1, quantile, probs = 0.9)

NI.sim.out.250.10 <- NI.sim.out[names(NI.sim.out.10) == "250"]
NI.sim.out.250.10 <- do.call(cbind, NI.sim.out.250.10)
NI.sim.out.250.10.mean <- rowMeans(NI.sim.out.250.10)
NI.sim.out.250.10.lcl  <- apply(NI.sim.out.250.10, 1, quantile, probs = 0.1)
NI.sim.out.250.10.ucl  <- apply(NI.sim.out.250.10, 1, quantile, probs = 0.9)

#15
NI.sim.out.150.15 <- NI.sim.out[names(NI.sim.out.15) == "150"]
NI.sim.out.150.15 <- do.call(cbind, NI.sim.out.150.15)
NI.sim.out.150.15.mean <- rowMeans(NI.sim.out.150.15)
NI.sim.out.150.15.lcl  <- apply(NI.sim.out.150.15, 1, quantile, probs = 0.1)
NI.sim.out.150.15.ucl  <- apply(NI.sim.out.150.15, 1, quantile, probs = 0.9)

NI.sim.out.200.15 <- NI.sim.out[names(NI.sim.out.15) == "200"]
NI.sim.out.200.15 <- do.call(cbind, NI.sim.out.200.15)
NI.sim.out.200.15.mean <- rowMeans(NI.sim.out.200.15)
NI.sim.out.200.15.lcl  <- apply(NI.sim.out.200.15, 1, quantile, probs = 0.1)
NI.sim.out.200.15.ucl  <- apply(NI.sim.out.200.15, 1, quantile, probs = 0.9)

NI.sim.out.250.15 <- NI.sim.out[names(NI.sim.out.15) == "250"]
NI.sim.out.250.15 <- do.call(cbind, NI.sim.out.250.15)
NI.sim.out.250.15.mean <- rowMeans(NI.sim.out.250.15)
NI.sim.out.250.15.lcl  <- apply(NI.sim.out.250.15, 1, quantile, probs = 0.1)
NI.sim.out.250.15.ucl  <- apply(NI.sim.out.250.15, 1, quantile, probs = 0.9)

#20
NI.sim.out.150.20 <- NI.sim.out[names(NI.sim.out.20) == "150"]
NI.sim.out.150.20 <- do.call(cbind, NI.sim.out.150.20)
NI.sim.out.150.20.mean <- rowMeans(NI.sim.out.150.20)
NI.sim.out.150.20.lcl  <- apply(NI.sim.out.150.20, 1, quantile, probs = 0.1)
NI.sim.out.150.20.ucl  <- apply(NI.sim.out.150.20, 1, quantile, probs = 0.9)

NI.sim.out.200.20 <- NI.sim.out[names(NI.sim.out.20) == "200"]
NI.sim.out.200.20 <- do.call(cbind, NI.sim.out.200.20)
NI.sim.out.200.20.mean <- rowMeans(NI.sim.out.200.20)
NI.sim.out.200.20.lcl  <- apply(NI.sim.out.200.20, 1, quantile, probs = 0.1)
NI.sim.out.200.20.ucl  <- apply(NI.sim.out.200.20, 1, quantile, probs = 0.9)

NI.sim.out.250.20 <- NI.sim.out[names(NI.sim.out.20) == "250"]
NI.sim.out.250.20 <- do.call(cbind, NI.sim.out.250.20)
NI.sim.out.250.20.mean <- rowMeans(NI.sim.out.250.20)
NI.sim.out.250.20.lcl  <- apply(NI.sim.out.250.20, 1, quantile, probs = 0.1)
NI.sim.out.250.20.ucl  <- apply(NI.sim.out.250.20, 1, quantile, probs = 0.9)

#25
NI.sim.out.150.25 <- NI.sim.out[names(NI.sim.out.25) == "150"]
NI.sim.out.150.25 <- do.call(cbind, NI.sim.out.150.25)
NI.sim.out.150.25.mean <- rowMeans(NI.sim.out.150.25)
NI.sim.out.150.25.lcl  <- apply(NI.sim.out.150.25, 1, quantile, probs = 0.1)
NI.sim.out.150.25.ucl  <- apply(NI.sim.out.150.25, 1, quantile, probs = 0.9)

NI.sim.out.200.25 <- NI.sim.out[names(NI.sim.out.25) == "200"]
NI.sim.out.200.25 <- do.call(cbind, NI.sim.out.200.25)
NI.sim.out.200.25.mean <- rowMeans(NI.sim.out.200.25)
NI.sim.out.200.25.lcl  <- apply(NI.sim.out.200.25, 1, quantile, probs = 0.1)
NI.sim.out.200.25.ucl  <- apply(NI.sim.out.200.25, 1, quantile, probs = 0.9)

NI.sim.out.250.25 <- NI.sim.out[names(NI.sim.out.25) == "250"]
NI.sim.out.250.25 <- do.call(cbind, NI.sim.out.250.25)
NI.sim.out.250.25.mean <- rowMeans(NI.sim.out.250.25)
NI.sim.out.250.25.lcl  <- apply(NI.sim.out.250.25, 1, quantile, probs = 0.1)
NI.sim.out.250.25.ucl  <- apply(NI.sim.out.250.25, 1, quantile, probs = 0.9)

#30
NI.sim.out.150.30 <- NI.sim.out[names(NI.sim.out.30) == "150"]
NI.sim.out.150.30 <- do.call(cbind, NI.sim.out.150.30)
NI.sim.out.150.30.mean <- rowMeans(NI.sim.out.150.30)
NI.sim.out.150.30.lcl  <- apply(NI.sim.out.150.30, 1, quantile, probs = 0.1)
NI.sim.out.150.30.ucl  <- apply(NI.sim.out.150.30, 1, quantile, probs = 0.9)

NI.sim.out.200.30 <- NI.sim.out[names(NI.sim.out.30) == "200"]
NI.sim.out.200.30 <- do.call(cbind, NI.sim.out.200.30)
NI.sim.out.200.30.mean <- rowMeans(NI.sim.out.200.30)
NI.sim.out.200.30.lcl  <- apply(NI.sim.out.200.30, 1, quantile, probs = 0.1)
NI.sim.out.200.30.ucl  <- apply(NI.sim.out.200.30, 1, quantile, probs = 0.9)

NI.sim.out.250.30 <- NI.sim.out[names(NI.sim.out.30) == "250"]
NI.sim.out.250.30 <- do.call(cbind, NI.sim.out.250.30)
NI.sim.out.250.30.mean <- rowMeans(NI.sim.out.250.30)
NI.sim.out.250.30.lcl  <- apply(NI.sim.out.250.30, 1, quantile, probs = 0.1)
NI.sim.out.250.30.ucl  <- apply(NI.sim.out.250.30, 1, quantile, probs = 0.9)

#35
NI.sim.out.150.35 <- NI.sim.out[names(NI.sim.out.35) == "150"]
NI.sim.out.150.35 <- do.call(cbind, NI.sim.out.150.35)
NI.sim.out.150.35.mean <- rowMeans(NI.sim.out.150.35)
NI.sim.out.150.35.lcl  <- apply(NI.sim.out.150.35, 1, quantile, probs = 0.1)
NI.sim.out.150.35.ucl  <- apply(NI.sim.out.150.35, 1, quantile, probs = 0.9)

NI.sim.out.200.35 <- NI.sim.out[names(NI.sim.out.35) == "200"]
NI.sim.out.200.35 <- do.call(cbind, NI.sim.out.200.35)
NI.sim.out.200.35.mean <- rowMeans(NI.sim.out.200.35)
NI.sim.out.200.35.lcl  <- apply(NI.sim.out.200.35, 1, quantile, probs = 0.1)
NI.sim.out.200.35.ucl  <- apply(NI.sim.out.200.35, 1, quantile, probs = 0.9)

NI.sim.out.250.35 <- NI.sim.out[names(NI.sim.out.35) == "250"]
NI.sim.out.250.35 <- do.call(cbind, NI.sim.out.250.35)
NI.sim.out.250.35.mean <- rowMeans(NI.sim.out.250.35)
NI.sim.out.250.35.lcl  <- apply(NI.sim.out.250.35, 1, quantile, probs = 0.1)
NI.sim.out.250.35.ucl  <- apply(NI.sim.out.250.35, 1, quantile, probs = 0.9)

#40
NI.sim.out.150.40 <- NI.sim.out[names(NI.sim.out.40) == "150"]
NI.sim.out.150.40 <- do.call(cbind, NI.sim.out.150.40)
NI.sim.out.150.40.mean <- rowMeans(NI.sim.out.150.40)
NI.sim.out.150.40.lcl  <- apply(NI.sim.out.150.40, 1, quantile, probs = 0.1)
NI.sim.out.150.40.ucl  <- apply(NI.sim.out.150.40, 1, quantile, probs = 0.9)

NI.sim.out.200.40 <- NI.sim.out[names(NI.sim.out.40) == "200"]
NI.sim.out.200.40 <- do.call(cbind, NI.sim.out.200.40)
NI.sim.out.200.40.mean <- rowMeans(NI.sim.out.200.40)
NI.sim.out.200.40.lcl  <- apply(NI.sim.out.200.40, 1, quantile, probs = 0.1)
NI.sim.out.200.40.ucl  <- apply(NI.sim.out.200.40, 1, quantile, probs = 0.9)

NI.sim.out.250.40 <- NI.sim.out[names(NI.sim.out.40) == "250"]
NI.sim.out.250.40 <- do.call(cbind, NI.sim.out.250.40)
NI.sim.out.250.40.mean <- rowMeans(NI.sim.out.250.40)
NI.sim.out.250.40.lcl  <- apply(NI.sim.out.250.40, 1, quantile, probs = 0.1)
NI.sim.out.250.40.ucl  <- apply(NI.sim.out.250.40, 1, quantile, probs = 0.9)

# group data
NI.sim.out.df.long <- data.frame(years = rep(years, times = 18),
                                 initN = rep(NI.pop.size, each = 31, times = 6),
                                 num_clutches = rep(c(10, 15, 20, 25, 30, 35),
                                                    each = 31, times = 3),
                                 mean = c(NI.sim.out.150.10.mean,
                                          NI.sim.out.200.10.mean,
                                          NI.sim.out.250.10.mean,
                                          NI.sim.out.150.15.mean,
                                          NI.sim.out.200.15.mean,
                                          NI.sim.out.250.15.mean,
                                          NI.sim.out.150.20.mean,
                                          NI.sim.out.200.20.mean,
                                          NI.sim.out.250.20.mean,
                                          NI.sim.out.150.25.mean,
                                          NI.sim.out.200.25.mean,
                                          NI.sim.out.250.25.mean,
                                          NI.sim.out.150.30.mean,
                                          NI.sim.out.200.30.mean,
                                          NI.sim.out.250.30.mean,
                                          NI.sim.out.150.35.mean,
                                          NI.sim.out.200.35.mean,
                                          NI.sim.out.250.35.mean),
                                          # NI.sim.out.150.40.mean,
                                          # NI.sim.out.200.40.mean,
                                          # NI.sim.out.250.40.mean),
                                 lcl = c(NI.sim.out.150.10.lcl,
                                         NI.sim.out.200.10.lcl,
                                         NI.sim.out.250.10.lcl,
                                         NI.sim.out.150.15.lcl,
                                         NI.sim.out.200.15.lcl,
                                         NI.sim.out.250.15.lcl,
                                         NI.sim.out.150.20.lcl,
                                         NI.sim.out.200.20.lcl,
                                         NI.sim.out.250.20.lcl,
                                         NI.sim.out.150.25.lcl,
                                         NI.sim.out.200.25.lcl,
                                         NI.sim.out.250.25.lcl,
                                         NI.sim.out.150.30.lcl,
                                         NI.sim.out.200.30.lcl,
                                         NI.sim.out.250.30.lcl,
                                         NI.sim.out.150.35.lcl,
                                         NI.sim.out.200.35.lcl,
                                         NI.sim.out.250.35.lcl),
                                         # NI.sim.out.150.40.lcl,
                                         # NI.sim.out.200.40.lcl,
                                         # NI.sim.out.250.40.lcl),
                                 ucl = c(NI.sim.out.150.10.ucl,
                                         NI.sim.out.200.10.ucl,
                                         NI.sim.out.250.10.ucl,
                                         NI.sim.out.150.15.ucl,
                                         NI.sim.out.200.15.ucl,
                                         NI.sim.out.250.15.ucl,
                                         NI.sim.out.150.20.ucl,
                                         NI.sim.out.200.20.ucl,
                                         NI.sim.out.250.20.ucl,
                                         NI.sim.out.150.25.ucl,
                                         NI.sim.out.200.25.ucl,
                                         NI.sim.out.250.25.ucl,
                                         NI.sim.out.150.30.ucl,
                                         NI.sim.out.200.30.ucl,
                                         NI.sim.out.250.30.ucl,
                                         NI.sim.out.150.35.ucl,
                                         NI.sim.out.200.35.ucl,
                                         NI.sim.out.250.35.ucl)
                                         # NI.sim.out.150.40.ucl,
                                         # NI.sim.out.200.40.ucl,
                                         # NI.sim.out.250.40.ucl)
)



write.table(NI.sim.out.df.long,
            "data/NI_A2_diffpops_clutches_1.csv",
            sep = ",",
            row.names = FALSE)


# Full sims plot
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
  xlab("Years") +
  ylab("Curlew population size (individuals)") +
  ylim(0, 900) +
  scale_x_continuous(breaks = seq(2022, 2052, 1)) +
  labs(colour = "Initial Starting Pop (breeding pairs)") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_bw() + 
  theme(axis.title = element_text(size = 32),
        axis.text = element_text(size = 16, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        legend.title = element_text(size = 20),
        legend.position = "top",
        legend.justification = "right",
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20))

NI.sim.out.plot
NI.sim.out.plot.ribbon



# Lower than x
bettermatch.100 <- which(NI.sim.out.100[31 ,] < 1,
                     arr.ind = TRUE)
NI.sim.out.100.match <- NI.sim.out.100[, bettermatch.100]
NI.sim.out.100.match.mean <- rowMeans(NI.sim.out.100.match)
NI.sim.out.100.match.lcl <- apply(NI.sim.out.100.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.100.match.ucl <- apply(NI.sim.out.100.match, 1,
                                   quantile, probs = 0.9)

bettermatch.150 <- which(NI.sim.out.150[31 ,] < 1,
                     arr.ind = TRUE)
NI.sim.out.150.match <- NI.sim.out.150[, bettermatch.150]
NI.sim.out.150.match.mean <- rowMeans(NI.sim.out.150.match)
NI.sim.out.150.match.lcl <- apply(NI.sim.out.150.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.150.match.ucl <- apply(NI.sim.out.150.match, 1,
                                   quantile, probs = 0.9)

bettermatch.200 <- which(NI.sim.out.200[31 ,] < 1,
                         arr.ind = TRUE)
NI.sim.out.200.match <- NI.sim.out.200[, bettermatch.200]
NI.sim.out.200.match.mean <- rowMeans(NI.sim.out.200.match)
NI.sim.out.200.match.lcl <- apply(NI.sim.out.200.match, 1,
                                  quantile, probs = 0.1)
NI.sim.out.200.match.ucl <- apply(NI.sim.out.200.match, 1,
                                  quantile, probs = 0.9)

bettermatch.250 <- which(NI.sim.out.250[31 ,] < 1,
                         arr.ind = TRUE)
NI.sim.out.250.match <- NI.sim.out.250[, bettermatch.250]
NI.sim.out.250.match.mean <- rowMeans(NI.sim.out.250.match)
NI.sim.out.250.match.lcl <- apply(NI.sim.out.250.match, 1,
                                  quantile, probs = 0.1)
NI.sim.out.250.match.ucl <- apply(NI.sim.out.250.match, 1,
                                  quantile, probs = 0.9)



NI.sim.out.df.match <- data.frame(years = rep(years, times = 4),
                                 initN = rep(NI.pop.size, each = 31),
                                 mean = c(NI.sim.out.100.match.mean,
                                          NI.sim.out.150.match.mean,
                                          NI.sim.out.200.match.mean,
                                          NI.sim.out.250.match.mean),
                                 lcl = c(NI.sim.out.100.match.lcl,
                                         NI.sim.out.150.match.lcl,
                                         NI.sim.out.200.match.lcl,
                                         NI.sim.out.250.match.lcl),
                                 ucl = c(NI.sim.out.100.match.ucl,
                                         NI.sim.out.150.match.ucl,
                                         NI.sim.out.200.match.ucl,
                                         NI.sim.out.250.match.ucl)
)


write.table(NI.sim.out.df.match,
            "data/NI_HS_A1_N100_hplr_means_match1_1.csv",
            sep = ",",
            row.names = FALSE)

# Plot less than 20
NI.sim.out.plot.match <- ggplot(NI.sim.out.df.match) +
  facet_wrap(~initN) +
  geom_line(aes(x = years,
                y = mean,
                colour = factor(initN)))

NI.sim.out.plot.ribbon.match <- NI.sim.out.plot.match +
  geom_ribbon(aes(x = years,
                  ymin = lcl,
                  ymax = ucl),
              alpha = 0.3) +
  xlab("Years") +
  ylab("Curlew population size (individuals)") +
  ylim(0, 600) +
  scale_x_continuous(breaks = seq(2022, 2052, 1)) +
  labs(colour = "Parameter Combinations") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_bw() + 
  theme(axis.title = element_text(size = 32),
        axis.text = element_text(size = 16, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        legend.title = element_text(size = 20),
        legend.position = "top",
        legend.justification = "right",
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20))

NI.sim.out.plot.match
NI.sim.out.plot.ribbon.match



#----- Calculate Extinction Rates -----

ext.out <- data.frame()


for (s in 1:length(sim.out)) {
  
  out <- data.frame()
  
  for (t in extinct.cutoff) {
    
    out$ext.thresh <- t
    out$outcome <- extinction.bysim(sim.out[s], t)
    
  }
  
}


sum(apply(sim.out[[1]], 1, function(t) min(t)<100))/ncol(sim.out[[1]])



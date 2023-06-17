######  Head-Starting Curlew Population Viability Analysis #####

# Head-Starting Scenario B - courtesy of Becs Lee
# Underlying productivity increases each year, 1st year survival increases
#   in the long-term
#   1. Scenario B2 & B3: no head-starting
#       Baseline model. Upon the advice of SO, changed model params to reduce
#       number of values for a) computational efficiency and b) irrelevance
#       to model output.
#       The model will be split into 6 different models:
#         a - low philopatry, high return
#         b - medium phil, high return
#         c - high phil, high return
#         d - low phil, low return
#         e - med phil, low return
#         f - high phil, low return
#         Low return is still a high prob at 0.7

# Increase in rates per year for how long? 
#     5 years? 8 years?
#     how much do rates increase?
# Would it make sense to use different matrices for each year with an increase?
#     Wherever fledgling or immature survival and productivity are
#     multiply by 1.1


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
NI.pop.size   <- seq(100, 150, 50)      # breeding pairs

ROI.pop.size  <- seq(100, 200, 25)


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

hs.pop        <- 0                      # baseline - 0 HS

capt.fail     <- 0.2                    # The number of HS birds that
# fail during rearing & can't be
# released

intro.years   <- c(1, 2, 3, 4, 5)       # years in which hs birds will be
# introduced

phil.p        <- 0.3                    # probability of first years return
# to release site

return.p      <- 0.70                   # probability of breeding birds 
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

extinct.cutoff  <- c(2, 6, 10, 20)      # below this number, pop is extinct

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
  S.jw*phil.p, (S.jhs*capt.fail)*phil.p,  0,  0,                      0,
  0,          0,      S.imm*return.p,         0,                      0,
  0,          0,      0,                      S.sub*return.p,     S.ad*return.p
)

# in second year - enough conservation measures to boost rates by 5 %
bird.matrix.2 <- expression(
  # f_w       f_hs      imm                     sub                     adult  
  0,          0,    ((prod.rates*1.05)*0.5)*0.005, ((prod.rates*1.05)*0.5)*0.1, (prod.rates*1.05)*0.5,
  0,          0,     0,                       0,                      0,
  (S.jw*1.05)*phil.p, ((S.jhs*1.05)*capt.fail)*phil.p,  0,  0,                      0,
  0,          0,      (S.imm*1.05)*return.p,         0,                      0,
  0,          0,      0,                      S.sub*return.p,     S.ad*return.p
)

bird.matrix.3 <- expression(
  # f_w       f_hs      imm                     sub                     adult  
  0,          0,    ((prod.rates*1.1)*0.5)*0.005, ((prod.rates*1.1)*0.5)*0.1, (prod.rates*1.1)*0.5,
  0,          0,     0,                       0,                      0,
  (S.jw*1.1)*phil.p, ((S.jhs*1.1)*capt.fail)*phil.p,  0,  0,                      0,
  0,          0,      (S.imm*1.1)*return.p,         0,                      0,
  0,          0,      0,                      S.sub*return.p,     S.ad*return.p
)

bird.matrix.4 <- expression(
  # f_w       f_hs      imm                     sub                     adult  
  0,          0,    (prod.rates*0.5)*0.005, (prod.rates*0.5)*0.1, prod.rates*0.5,
  0,          0,     0,                       0,                      0,
  S.jw*phil.p, (S.jhs*capt.fail)*phil.p,  0,  0,                      0,
  0,          0,      S.imm*return.p,         0,                      0,
  0,          0,      0,                      S.sub*return.p,     S.ad*return.p
)

bird.matrix.5 <- expression(
  # f_w       f_hs      imm                     sub                     adult  
  0,          0,    (prod.rates*0.5)*0.005, (prod.rates*0.5)*0.1, prod.rates*0.5,
  0,          0,     0,                       0,                      0,
  S.jw*phil.p, (S.jhs*capt.fail)*phil.p,  0,  0,                      0,
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
                     "years.HS")

write.table(simul.in,
            "data/NI_HS_scenario_A1_base_N100_hplr_parameters.csv",
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
    Pop.Array[1, rep]   <- (init.N*2) + hs.numbers
    
    
    # Loop through the years
    for (y in 2:(n.years + 1)) {
      
      # create matrix with vital rates
      # Check if ruined breeding season - varies every sim
      F.Year <- ifelse(rbinom(1, 1, breed.fail) == 1,
                       F.loop * 0.5,
                       F.loop)
      
      # Create list of vital rates
      bird.vr <- list(prod.rates  = F.Year,
                      S.ad        = S.ad.loop,
                      S.jw        = S.jw.loop,
                      S.jhs       = S.jhs.loop,
                      S.imm       = S.imm.loop,
                      S.sub       = S.sub.loop,
                      phil.p      = philopatry,
                      return.p    = sitefidelity,
                      capt.fail   = capt.fail
      )
      
      # Create leslie matrix
      if (y == 1) {
        a <- matrix(sapply(bird.matrix,
                           eval,
                           bird.vr,
                           NULL),
                    nrow = sqrt(length(bird.matrix)),
                    byrow = TRUE)
      } else if (y == 2) {
        a <- matrix(sapply(bird.matrix.2,
                           eval,
                           bird.vr,
                           NULL),
                    nrow = sqrt(length(bird.matrix)),
                    byrow = TRUE)
      } else if (y == 3) {
        a <- matrix(sapply(bird.matrix.3,
                           eval,
                           bird.vr,
                           NULL),
                    nrow = sqrt(length(bird.matrix)),
                    byrow = TRUE)
      } else if (y == 4) {
        a <- matrix(sapply(bird.matrix.4,
                           eval,
                           bird.vr,
                           NULL),
                    nrow = sqrt(length(bird.matrix)),
                    byrow = TRUE)
      } else if (y == 5) {
        a <- matrix(sapply(bird.matrix.5,
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
      if (y == 2) {
        pop.size  <- c((F.loop * 0.5 * (init.N * 2)), # Jws
                       0,  # Jhs
                       ((hs.numbers * (1-capt.fail)) * S.jhs.loop * philopatry), # Imm
                       0, # sub
                       ((init.N * 2) * S.ad.loop * sitefidelity)) # adu
      } else if (y %in% intro.years) {
        pop.size <- c(projections$stage.vectors[1, y - 1],
                      (hs.numbers * (1 - capt.fail)),
                      projections$stage.vectors[3, y - 1],
                      projections$stage.vectors[4, y - 1],
                      projections$stage.vectors[5, y - 1])
      } else {
        pop.size <- c(projections$stage.vectors[, y - 1])
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
      if (y < hs.years) next.year <- next.year + hs.numbers
      
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
  
  NI.sim.out <- foreach(s = c(1:dim(simul.in)[1]),
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
saveRDS(NI.sim.out,
        "data/NI_HS_base_A1_N100_sims_hplr_1.RDS")



#----- Population plots -----
NI.sim.out <- readRDS("data/NI_HS_base_A1_sims_mplr_1.RDS")
names.NI.sim.out <- rep(as.character(c(NI.pop.size)), (22400/2))

names(NI.sim.out) <- names.NI.sim.out
NI.sim.out <- sim.out

years <- seq(2022, 2052, 1)

NI.sim.out.100 <- NI.sim.out[names(NI.sim.out) == "100"]
NI.sim.out.100 <- do.call(cbind, NI.sim.out.100)
NI.sim.out.100.mean <- rowMeans(NI.sim.out.100)
NI.sim.out.100.lcl  <- apply(NI.sim.out.100, 1, quantile, probs = 0.1)
NI.sim.out.100.ucl  <- apply(NI.sim.out.100, 1, quantile, probs = 0.9)

NI.sim.out.150 <- NI.sim.out[names(NI.sim.out) == "150"]
NI.sim.out.150 <- do.call(cbind, NI.sim.out.150)
NI.sim.out.150.mean <- rowMeans(NI.sim.out.150)
NI.sim.out.150.lcl  <- apply(NI.sim.out.150, 1, quantile, probs = 0.1)
NI.sim.out.150.ucl  <- apply(NI.sim.out.150, 1, quantile, probs = 0.9)

NI.sim.out.200 <- NI.sim.out[names(NI.sim.out) == "200"]
NI.sim.out.200 <- do.call(cbind, NI.sim.out.200)
NI.sim.out.200.mean <- rowMeans(NI.sim.out.200)
NI.sim.out.200.lcl  <- apply(NI.sim.out.200, 1, quantile, probs = 0.1)
NI.sim.out.200.ucl  <- apply(NI.sim.out.200, 1, quantile, probs = 0.9)

NI.sim.out.225 <- NI.sim.out[names(NI.sim.out) == "225"]
NI.sim.out.225 <- do.call(cbind, NI.sim.out.225)
NI.sim.out.225.mean <- rowMeans(NI.sim.out.225)
NI.sim.out.225.lcl  <- apply(NI.sim.out.225, 1, quantile, probs = 0.1)
NI.sim.out.225.ucl  <- apply(NI.sim.out.225, 1, quantile, probs = 0.9)

NI.sim.out.250 <- NI.sim.out[names(NI.sim.out) == "250"]
NI.sim.out.250 <- do.call(cbind, NI.sim.out.250)
NI.sim.out.250.mean <- rowMeans(NI.sim.out.250)
NI.sim.out.250.lcl  <- apply(NI.sim.out.250, 1, quantile, probs = 0.1)
NI.sim.out.250.ucl  <- apply(NI.sim.out.250, 1, quantile, probs = 0.9)

NI.sim.out.275 <- NI.sim.out[names(NI.sim.out) == "275"]
NI.sim.out.275 <- do.call(cbind, NI.sim.out.275)
NI.sim.out.275.mean <- rowMeans(NI.sim.out.275)
NI.sim.out.275.lcl  <- apply(NI.sim.out.275, 1, quantile, probs = 0.1)
NI.sim.out.275.ucl  <- apply(NI.sim.out.275, 1, quantile, probs = 0.9)

NI.sim.out.300 <- NI.sim.out[names(NI.sim.out) == "300"]
NI.sim.out.300 <- do.call(cbind, NI.sim.out.300)
NI.sim.out.300.mean <- rowMeans(NI.sim.out.300)
NI.sim.out.300.lcl  <- apply(NI.sim.out.300, 1, quantile, probs = 0.1)
NI.sim.out.300.ucl  <- apply(NI.sim.out.300, 1, quantile, probs = 0.9)


NI.sim.out.df.long <- data.frame(years = rep(years, times = 2),
                                 initN = rep(NI.pop.size, each = 31),
                                 mean = c(NI.sim.out.100.mean,
                                          NI.sim.out.150.mean),
                                 lcl = c(NI.sim.out.100.lcl,
                                         NI.sim.out.150.lcl),
                                 ucl = c(NI.sim.out.100.ucl,
                                         NI.sim.out.150.ucl)
)


write.table(NI.sim.out.df.long,
            "data/NI_HS_base_A1_N100_means_hplr_fullsims_1.csv",
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
  ylim(0, 600) +
  theme_bw() + 
  theme(axis.title = element_text(size = 32),
        axis.text = element_text(size = 16, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        legend.position = "none")

NI.sim.out.plot
NI.sim.out.plot.ribbon

# Lower than x
bettermatch.100 <- which(NI.sim.out.100[31 ,] < 5,
                         arr.ind = TRUE)
NI.sim.out.100.match <- NI.sim.out.100[, bettermatch.100]
NI.sim.out.100.match.mean <- rowMeans(NI.sim.out.100.match)
NI.sim.out.100.match.lcl <- apply(NI.sim.out.100.match, 1,
                                  quantile, probs = 0.1)
NI.sim.out.100.match.ucl <- apply(NI.sim.out.100.match, 1,
                                  quantile, probs = 0.9)

bettermatch.150 <- which(NI.sim.out.150[31 ,] < 5,
                         arr.ind = TRUE)
NI.sim.out.150.match <- NI.sim.out.150[, bettermatch.150]
NI.sim.out.150.match.mean <- rowMeans(NI.sim.out.150.match)
NI.sim.out.150.match.lcl <- apply(NI.sim.out.150.match, 1,
                                  quantile, probs = 0.1)
NI.sim.out.150.match.ucl <- apply(NI.sim.out.150.match, 1,
                                  quantile, probs = 0.9)



NI.sim.out.df.match <- data.frame(years = rep(years, times = 2),
                                  initN = rep(NI.pop.size, each = 31),
                                  mean = c(NI.sim.out.100.match.mean,
                                           NI.sim.out.150.match.mean),
                                  lcl = c(NI.sim.out.100.match.lcl,
                                          NI.sim.out.150.match.lcl),
                                  ucl = c(NI.sim.out.100.match.ucl,
                                          NI.sim.out.150.match.ucl)
)


write.table(NI.sim.out.df.match,
            "data/NI_HS_base_A1_N100_hplr_means_match5_1.csv",
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
  ylim(0, 500) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 32),
    axis.text = element_text(size = 16, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(),
    legend.position = "none")

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



######  Head-Starting Curlew Population Viability Analysis #####

# Head-Starting Scenario A - courtesy of Becs Lee
# Underlying productivity does not change, low productivity continues
#   in the long-term
#   1. Scenario Individual Sites:
#       - testing the population response to removing clutches from nests
#         for 5 consecutive years - limited relay probability
#       - no head-starting
#       - need to run simplified models for bands of small pops
#       -   1-4, 5-9, 10-14 - just do 2, 4, 10, 14
#       = glenwherry example of 58 breeding pairs
# Also model the number of clutches required to boost the population
#  by 20 breeding pairs within 10 years

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
#NI.pop.size   <- seq(100, 250, 50)      # breeding pairs
site.pop      <- 4                     # number of breed pairs on site
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

hs.clutches   <- 0
hs.pop        <- (4*hs.clutches)        # baseline - 0 HS
# capt fail = 20% - start with 10 clutches then increase

capt.fail     <- 0.2                    # The number of HS birds that
# fail during rearing & can't be
# released

intro.years   <- 5           # years in which hs birds will be
# introduced

phil.p        <- 0.4                    # probability of first years return
# to release site
# philopatry is still possible in wild
# bred birds, but not HS

return.p      <- 0.99                   # probability of breeding birds 
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

breed.fail    <- 0.4                    # prob the breeding season is a 
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
  0,          0,    (prod.rates*0.5)*0.005, (prod.rates*0.5)*0.1, (prod.rates*1.2)*0.5,
  0,          0,     0,                       0,                      0,
  S.jw*phil.p, (S.jhs*capt.fail)*phil.p,  0,  0,                      0,
  0,          0,      S.imm*return.p,         0,                      0,
  0,          0,      0,                      S.sub*return.p,     S.ad*return.p
)

bird.matrix.5 <- expression(
  # f_w       f_hs      imm                     sub                     adult  
  0,          0,    (prod.rates*0.5)*0.005, (prod.rates*0.5)*0.1, (prod.rates*1.5)*0.5,
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
  (S.jw*phil.p)*0, (S.jhs*capt.fail)*(phil.p*0),  0,  0,                      0,
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
            "data/curlew_HS_NI_simple_params_B1.csv",
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

testPVA <- PVA(
  n.reps,
  n.years,
  init.N = simul.in.simple[s, 1],
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

testPVA.means <- as.data.frame(rowMeans(testPVA))
testPVA.lcl  <- apply(testPVA, 1, quantile, probs = 0.1)
testPVA.ucl  <- apply(testPVA, 1, quantile, probs = 0.9)
# What to name each element? intro bird numbers?


PVA_min <- testPVA
PVA_med <- testPVA
PVA_max <- testPVA

years <- seq(2023, 2053, 1)


PVA_ROI_A1 <- as.data.frame(cbind(
  years,
  PVA_min,
  PVA_med,
  PVA_max
))

PVA_ROI_A1_means <- as.data.frame(cbind(
  years,
  testPVA.means,
  testPVA.lcl,
  testPVA.ucl
))

PVA_ROI_A2 <- as.data.frame(cbind(
  years,
  PVA_min,
  PVA_med,
  PVA_max
))

PVA_ROI_A3 <- as.data.frame(cbind(
  years,
  PVA_min,
  PVA_med,
  PVA_max
))

PVA_ROI_A4 <- as.data.frame(cbind(
  years,
  PVA_min,
  PVA_med,
  PVA_max
))



ROI_A_pva <- as.tibble(rbind(
  PVA_ROI_A1,
  PVA_ROI_A2,
  PVA_ROI_A3,
  PVA_ROI_A4
))

ROI_A_means <- as.tibble(cbind(
  rowMeans(PVA_ROI_A1[, -1])/2,
  rowMeans(PVA_ROI_A2[, -1])/2,
  rowMeans(PVA_ROI_A3[, -1])/2,
  rowMeans(PVA_ROI_A4[, -1])/2
))

ROI_A_means$years <- years
names(ROI_A_means) <- c("A1", "A2", "A3", "A4", "years")




num_clutches <- "num_clutches"
scenario <- "scenario"
pop_size <- "population size (breeding pairs)"
gathercols <- c("A1", "A2", "A3", "A4")

ROI_A <- gather(ROI_A_means,
               scenario, pop_size, gathercols)


write.table(ROI_A,
            "data/HS_ROI_A_means_simple_1.csv",
            sep = ",",
            row.names = FALSE)

ROI_A__max_plot <- ggplot(ROI_A_max) +
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            size = 1) +
  # geom_smooth(aes(x = years, y = pop_size,
  #                 group = scenario,
  #                 colour = scenario),
  #             size = 1) +
  scale_colour_manual(values = foxPal) +
  xlab("Years") +
  ylab("Curlew population size (breeding pairs)") +
  ylim(0, 300) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(2023, 2054),
                     breaks = seq(2023, 2053, 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), 
                     limits = c(0, 300)) +
  labs(colour = "Scenario") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme_bw() + 
  theme(axis.title = element_text(size = 32,
                                  vjust = -0.1),
        axis.text = element_text(size = 20, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        legend.title = element_text(size = 20),
        legend.position = c(0.9, 0.9),
        legend.justification = "right",
        legend.text = element_text(size = "20"),
        strip.text.x = element_text(size = 20))

PVA_min <- testPVA
PVA_med <- testPVA
PVA_max <- testPVA

PVA_ROI_B1 <- as.data.frame(cbind(
  years,
  PVA_min,
  PVA_med,
  PVA_max
))


PVA_ROI_B2 <- as.data.frame(cbind(
  years,
  PVA_min,
  PVA_med,
  PVA_max
))

PVA_ROI_B3 <- as.data.frame(cbind(
  years,
  PVA_min,
  PVA_med,
  PVA_max
))

PVA_ROI_B4 <- as.data.frame(cbind(
  years,
  PVA_min,
  PVA_med,
  PVA_max
))



ROI_B_pva <- as.tibble(rbind(
  PVA_ROI_B1,
  PVA_ROI_B2,
  PVA_ROI_B3,
  PVA_ROI_B4
))

ROI_B_means <- as.tibble(cbind(
  rowMeans(PVA_ROI_B1[, -1])/2,
  rowMeans(PVA_ROI_B2[, -1])/2,
  rowMeans(PVA_ROI_B3[, -1])/2,
  rowMeans(PVA_ROI_B4[, -1])/2
))

ROI_B_means$years <- years
names(ROI_B_means) <- c("B1", "B2", "B3", "B4", "years")


scenario <- "scenario"
pop_size <- "population size (breeding pairs)"
gathercols <- c("B1", "B2", "B3", "B4")

ROI_B <- gather(ROI_B_means,
               scenario, pop_size, gathercols)




write.table(ROI_B,
            "data/HS_ROI_B_means_simple_1.csv",
            sep = ",",
            row.names = FALSE)

ROI_B_plot <- ggplot(ROI_B) +
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario)) +
  xlab("Years") +
  ylab("Curlew population size (breeding pairs)") +
  ylim(0, 300) +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(2023, 2054),
                     breaks = seq(2023, 2053, 1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 300)) +
  labs(colour = "Scenario") +
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
        legend.text = element_text(size = "20"),
        strip.text.x = element_text(size = 20))

glen20bp_facet <- ggplot(glen20bp) +
  facet_wrap(~num_clutches) +
  geom_line(aes(x = years,
                y = pop_size,
                group = values,
                colour = values)) +
  stat_summary(aes(x = years,
                   y = pop_size,
                   group = initN),
               geom = "smooth",
               fun = "mean",
               colour = "black") +
  xlab("Years") +
  ylab("Curlew population size (individuals)") +
  ylim(0, 200) +
  scale_x_continuous(breaks = seq(2022, 2052, 1)) +
  labs(colour = "Parameters (min, mean, max)") +
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
        legend.text = element_text(size = "20"),
        strip.text.x = element_text(size = 20))






#######################

siteRemPVA <- as_tibble(rbind(
  PVA_lower,
  PVA_middle,
  PVA_upper,
  PVA_20
))


siteRemPVA <- rbind(PVA_glen,
                    PVA_lower,
                    PVA_middle,
                    PVA_upper,
                    PVA_20)
names(siteRemPVA) <- c("years", "low values", "med values", "high values")
values <- "values"
pop_size <- "population size (breeding pairs)"
gathercols <- c("low values", "med values", "high values")

siteRemPVA_long <- gather(siteRemPVA,
                          values, pop_size, gathercols)
initN <- rep(c(4, 8, 12, 20), each = 31)

glen_long <- gather(PVA_glen,
                    values,
                    pop_size,
                    gathercols)
glen_long$initN <- rep(c(58), each = 93)

siteRemPVA_long$initN <- rep(initN, times = 3) 

write.table(siteRemPVA_long,
            "data/HS_site_Rem_PVA_bands_lplr_1.csv",
            sep = ",",
            row.names = FALSE)

# Plot
siteRem_plot <- ggplot(siteRemPVA_long) +
  geom_line(aes(x = years,
                y = pop_size,
                group = values,
                colour = values)) +
  xlab("Years") +
  ylab("Curlew Pop Size (individs)") +
  ylim(0, 300) +
  scale_x_continuous(breaks = seq(2022, 2052, 1)) +
  theme_bw() +
  theme(axis.title = element_text(size = 32),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(angle = 4, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        legend.position = "none")

siteRem_plot_facet <- ggplot(glen_long) +
  facet_wrap(~initN) +
  geom_line(aes(x = years,
                y = pop_size,
                group = values,
                colour = values)) +
  stat_summary(aes(x = years,
                   y = pop_size,
                   group = initN),
               geom = "smooth",
               fun = "mean",
               colour = "black") +
  xlab("Years") +
  ylab("Curlew population size (individuals)") +
  ylim(0, 200) +
  scale_x_continuous(breaks = seq(2022, 2052, 1)) +
  labs(colour = "Parameters (min, mean, max)") +
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
        legend.text = element_text(size = "20"),
        strip.text.x = element_text(size = 20))
siteRem_plot_facet




# Write to save file
saveRDS(ROI.sim.out,
        "data/ROI_HS_base_A4_sims_lplr_1.RDS")



#----- Population plots -----
ROI.sim.out <- readRDS("data/ROI_HS_A4_sims_lplr_1.RDS")
names.ROI.sim.out <- rep(as.character(c(ROI.pop.size)), (33600/3))

names(ROI.sim.out) <- names.ROI.sim.out
ROI.sim.out <- sim.out

years <- seq(2022, 2052, 1)

ROI.sim.out.100 <- ROI.sim.out[names(ROI.sim.out) == "100"]
ROI.sim.out.100 <- do.call(cbind, ROI.sim.out.100)
ROI.sim.out.100.mean <- rowMeans(ROI.sim.out.100)
ROI.sim.out.100.lcl  <- apply(ROI.sim.out.100, 1, quantile, probs = 0.1)
ROI.sim.out.100.ucl  <- apply(ROI.sim.out.100, 1, quantile, probs = 0.9)

ROI.sim.out.110 <- ROI.sim.out[names(ROI.sim.out) == "110"]
ROI.sim.out.110 <- do.call(cbind, ROI.sim.out.110)
ROI.sim.out.110.mean <- rowMeans(ROI.sim.out.110)
ROI.sim.out.110.lcl  <- apply(ROI.sim.out.110, 1, quantile, probs = 0.1)
ROI.sim.out.110.ucl  <- apply(ROI.sim.out.110, 1, quantile, probs = 0.9)

ROI.sim.out.120 <- ROI.sim.out[names(ROI.sim.out) == "120"]
ROI.sim.out.120 <- do.call(cbind, ROI.sim.out.120)
ROI.sim.out.120.mean <- rowMeans(ROI.sim.out.120)
ROI.sim.out.120.lcl  <- apply(ROI.sim.out.120, 1, quantile, probs = 0.1)
ROI.sim.out.120.ucl  <- apply(ROI.sim.out.120, 1, quantile, probs = 0.9)


ROI.sim.out.df.long <- data.frame(years = rep(years, times = 3),
                                  initN = rep(ROI.pop.size, each = 31),
                                  mean = c(ROI.sim.out.100.mean,
                                           ROI.sim.out.110.mean,
                                           ROI.sim.out.120.mean),
                                  lcl = c(ROI.sim.out.100.lcl,
                                          ROI.sim.out.110.lcl,
                                          ROI.sim.out.120.lcl),
                                  ucl = c(ROI.sim.out.100.ucl,
                                          ROI.sim.out.110.ucl,
                                          ROI.sim.out.120.ucl)
)


write.table(ROI.sim.out.df.long,
            "data/ROI_HS_base_A4_means_lplr_fullsims_1.csv",
            sep = ",",
            row.names = FALSE)


# Full sims plot
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
  xlab("Years") +
  ylab("Curlew population size (individuals)") +
  ylim(0, 500) +
  scale_x_continuous(breaks = seq(2022, 2052, 1)) +
  labs(colour = "Initial starting population (breeding pairs)") +
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
        legend.text = element_text(size = "20"),
        strip.text.x = element_text(size = 20))

ROI.sim.out.plot
ROI.sim.out.plot.ribbon

# Lower than x
bettermatch.100 <- which(ROI.sim.out.100[31 ,] < 5,
                         arr.ind = TRUE)
ROI.sim.out.100.match <- ROI.sim.out.100[, bettermatch.100]
ROI.sim.out.100.match.mean <- rowMeans(ROI.sim.out.100.match)
ROI.sim.out.100.match.lcl <- apply(ROI.sim.out.100.match, 1,
                                   quantile, probs = 0.1)
ROI.sim.out.100.match.ucl <- apply(ROI.sim.out.100.match, 1,
                                   quantile, probs = 0.9)

bettermatch.110 <- which(ROI.sim.out.110[31 ,] < 5,
                         arr.ind = TRUE)
ROI.sim.out.110.match <- ROI.sim.out.110[, bettermatch.110]
ROI.sim.out.110.match.mean <- rowMeans(ROI.sim.out.110.match)
ROI.sim.out.110.match.lcl <- apply(ROI.sim.out.110.match, 1,
                                   quantile, probs = 0.1)
ROI.sim.out.110.match.ucl <- apply(ROI.sim.out.110.match, 1,
                                   quantile, probs = 0.9)

bettermatch.120 <- which(ROI.sim.out.120[31 ,] < 5,
                         arr.ind = TRUE)
ROI.sim.out.120.match <- ROI.sim.out.120[, bettermatch.120]
ROI.sim.out.120.match.mean <- rowMeans(ROI.sim.out.120.match)
ROI.sim.out.120.match.lcl <- apply(ROI.sim.out.120.match, 1,
                                   quantile, probs = 0.1)
ROI.sim.out.120.match.ucl <- apply(ROI.sim.out.120.match, 1,
                                   quantile, probs = 0.9)



ROI.sim.out.df.match <- data.frame(years = rep(years, times = 2),
                                   initN = rep(ROI.pop.size, each = 31),
                                   mean = c(ROI.sim.out.100.match.mean,
                                            ROI.sim.out.110.match.mean,
                                            ROI.sim.out.120.mattch.mean),
                                   lcl = c(ROI.sim.out.100.match.lcl,
                                           ROI.sim.out.110.match.lcl,
                                           ROI.sim.out.120.match.lcl),
                                   ucl = c(ROI.sim.out.100.match.ucl,
                                           ROI.sim.out.110.match.ucl,
                                           ROI.sim.out.120.match.ucl)
)


write.table(ROI.sim.out.df.match,
            "data/ROI_HS_base_A4_hplr_means_match1_1.csv",
            sep = ",",
            row.names = FALSE)

# Plot less than 20
ROI.sim.out.plot.match <- ggplot(ROI.sim.out.df.match) +
  facet_wrap(~initN) +
  geom_line(aes(x = years,
                y = mean,
                colour = factor(initN)))

ROI.sim.out.plot.ribbon.match <- ROI.sim.out.plot.match +
  geom_ribbon(aes(x = years,
                  ymin = lcl,
                  ymax = ucl),
              alpha = 0.3) +
  xlab("Years") +
  ylab("Curlew population size (individuals)") +
  ylim(0, 500) +
  scale_x_continuous(breaks = seq(2022, 2052, 1)) +
  labs(colour = "Initial starting population (breeding pairs)") +
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
        legend.text = element_text(size = "20"),
        strip.text.x = element_text(size = 20))

ROI.sim.out.plot.match
ROI.sim.out.plot.ribbon.match



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






#------- Write Files to Disk  ------
write.table(ROI_A_means,
            "data/curlew_HS_ROI_A_means_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_ROI_A1,
            "data/curlew_HS_ROI_A1_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_ROI_A2,
            "data/curlew_HS_ROI_A2_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_ROI_A3,
            "data/curlew_HS_ROI_A3_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_ROI_A4,
            "data/curlew_HS_ROI_A4_simple_1.csv",
            sep = ",",
            row.names = FALSE)



# NI
write.table(PVA_NI_A1,
            "data/curlew_HS_NI_A1_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_NI_A2,
            "data/curlew_HS_NI_A2_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_NI_A3,
            "data/curlew_HS_NI_A3_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_NI_A4,
            "data/curlew_HS_NI_A4_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_NI_B1,
            "data/curlew_HS_NI_B1_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_NI_B2,
            "data/curlew_HS_NI_B2_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_NI_B3,
            "data/curlew_HS_NI_B3_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_NI_B4,
            "data/curlew_HS_NI_B4_simple_1.csv",
            sep = ",",
            row.names = FALSE)


write.table(NI_A_means,
            "data/curlew_HS_NI_A_means_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(NI_B_means,
            "data/curlew_HS_NI_B_means_simple_1.csv",
            sep = ",",
            row.names = FALSE)


# Glen & LLE
write.table(PVA_glen,
            "data/curlew_HS_glen_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_glen_rem,
            "data/curlew_HS_glen_rem_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_glen_10,
            "data/curlew_HS_glen_10_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_glen_15,
            "data/curlew_HS_glen_15_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_glen_20,
            "data/curlew_HS_glen_20_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_glen_25,
            "data/curlew_HS_glen_25_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_glen_30,
            "data/curlew_HS_glen_30_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_glen_35,
            "data/curlew_HS_glen_35_simple_1.csv",
            sep = ",",
            row.names = FALSE)


# LLE
write.table(PVA_lle,
            "data/curlew_HS_lle_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_lle_rem,
            "data/curlew_HS_lle_rem_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_lle_10,
            "data/curlew_HS_lle_10_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_lle_15,
            "data/curlew_HS_lle_15_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_lle_20,
            "data/curlew_HS_lle_20_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_lle_25,
            "data/curlew_HS_lle_25_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_lle_30,
            "data/curlew_HS_lle_30_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_lle_35,
            "data/curlew_HS_lle_35_simple_1.csv",
            sep = ",",
            row.names = FALSE)


# PVA Starting Pops
write.table(PVA_2,
            "data/curlew_HS_PVA2_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_2_rem,
            "data/curlew_HS_PVA2_rem_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_4,
            "data/curlew_HS_PVA4_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_4_rem,
            "data/curlew_HS_PVA4_rem_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_8,
            "data/curlew_HS_PVA8_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_8_rem,
            "data/curlew_HS_PVA8_rem_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_12,
            "data/curlew_HS_PVA12_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_12_rem,
            "data/curlew_HS_PVA12_rem_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_20,
            "data/curlew_HS_PVA20_simple_1.csv",
            sep = ",",
            row.names = FALSE)

write.table(PVA_20_rem,
            "data/curlew_HS_PVA20_rem_simple_1.csv",
            sep = ",",
            row.names = FALSE)



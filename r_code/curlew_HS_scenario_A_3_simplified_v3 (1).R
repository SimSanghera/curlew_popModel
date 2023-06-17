######  Head-Starting Curlew Population Viability Analysis #####

# Head-Starting Scenario A - courtesy of Becs Lee
# Underlying productivity does not change, low productivity continues
#   in the long-term
#   3. Scenario A3: maximum head-starting possible
#       This is based on the avaialblity and access to eggs within the local
#         area. From the NI table as it stands right now, this estimate is 
#         around 25 clutches per year.
#         With the hope that curlew re-lay.

# 21st Feb Update
# As an addition to the original population model, I am going to run 3
# separate models using the lowest, mean and highest values for each
# parameter to create a simple visual output similar to Becs' Excel model.
# Hopefully this will a) allow people to see something they can relate to,
# b) enable quicker changing of the headstarting birds parameter to determine
# the number of HS birds needed to stabilise the population.
# In the original model, this is difficult to see as we have so much
# uncertainty built into the parameters and model structure, obtaining a simple
# population time series plot is very difficult.
# The original model will be beneficial for identifying the probability of the 
# population going extinct under different HS scenarios and egg numbers.


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
NI.pop.size   <- c(100, 150, 200)      # breeding pairs
                                       # 07/03/22 Change init N to 100 and 150
site.pop      <- 12


# Survival Rates & Curlew Related Parameters
S.ad          <- c(0.78, 0.85, 0.92)  # Cook et al., 2021
# Roodbergen et al., 2021

S.jw          <- c(0.278, 0.326, 0.378)  # Wild bird survival
# Cook et al., 2021

S.jhs         <- c(0.278, 0.326, 0.378)  # HS bird survival
# Assume same as wild for now
# Geoff Hilton (pers comms)

S.sub         <- c(0.78, 0.85,  0.92)  # subadult survival

S.imm         <- c(0.40, 0.55, 0.70)  # survival of immature (age fledge - 1)

prod.rates    <- c(0.14, 0.25, 0.52)  # Grant et al., 1999
# Harris et al., 2019
# Zielonka et al., 2019
# Irish Breeding Bird Report

hs.clutches   <- 4                   # number of clutches - assume 4 eggs
hs.pop        <- (4*hs.clutches)        # maximum eggs feasible

capt.fail     <- 0.2                    # The number of HS birds that
# fail during rearing & can't be
#  released

intro.years   <- 5                       # years in which hs birds will be
# introduced

phil.p        <- c(0.3, 0.5, 0.7)       # probability of first years return
# to release site

return.p      <- c(0.3, 0.6, 0.90)      # probability of breeding birds 
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

breed.fail    <- 0.5                    # prob the breeding season is a 
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
  S.jw*phil.p, (S.jhs*(1-capt.fail)),  0,  0,                      0,
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
            "data/HS_scenario_A3_sites_simplified_parameters_1.csv",
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
        pop.size  <- c((F.loop * 0.5 * (init.N * 2)),
                       (hs.numbers * (1-capt.fail)),
                       ((hs.numbers * (1-capt.fail)) * S.jhs.loop * philopatry),
                       0,
                       ((init.N * 2) * S.ad.loop * sitefidelity))
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
                                    iterations  = n.years)
      
      # Stochastic Ricker Model
      R.max <- projections$lambda
      
      # Calculate abundance
      next.year   <- max(0,
                         trunc(Ricker(Pop.Array[y - 1, rep])))
      
      # Catastrophe
      if(runif(1) < catastrophe.p) next.year <- next.year * catastrophe.i
     
      # cat.p <- runif(1)
      # 
      # next.year <- ifelse(cat.p < catastrophe.p && y <= hs.years, 
      #                     next.year <- (next.year * catastrophe.i) + hs.numbers,
      #                     ifelse(cat.p < catastrophe.p && y > hs.years, 
      #                            next.year <- next.year * catastrophe.i,
      #                            next.year <- next.year))
      
      # Put in array
      Pop.Array[y, rep] <- next.year
      
    } # end year
    
  }   # end rep
  
  return(Pop.Array)
  
}


#---- Test PVA -----
testPVA <- PVA(
  n.reps,
  n.years,
  init.N = simul.in[s, 1],
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



years <- seq(2022, 2052, 1)

PVA_min <- testPVA
PVA_med <- testPVA
PVA_max <- testPVA

PVA_glen <- as.data.frame(cbind(
  years,
  PVA_min,
  PVA_med,
  PVA_max
))


PVA_lower_4clutch <- as.data.frame(cbind(
  years,
  PVA_min,
  PVA_med,
  PVA_max
))

PVA_lower_6clutch <- as.data.frame(cbind(
  years,
  PVA_min,
  PVA_med,
  PVA_max
))

PVA_lower_10clutch <- as.data.frame(cbind(
  years,
  PVA_min,
  PVA_med,
  PVA_max
))

PVA_lower_allclutch <- rbind(PVA_lower_4clutch,
                             PVA_lower_6clutch,
                             PVA_lower_10clutch)
PVA_lower_allclutch$clutchHS <- rep(c(4, 6, 10), each = 31)

names(PVA_lower_allclutch) <- c("years",
                             "low values",
                             "med values",
                             "high values",
                             "clutches HS")


values <- "values"
pop_size <- "population size (individs)"
gathercols <- c("low values", "med values", "high values")

PVA_low_all_long <- gather(PVA_lower_allclutch,
                          values, pop_size, gathercols)
initN <- rep(c(4, 8, 12), each = 31)

PVA_low_all_long$initN <- rep(initN, times = 3) 


PVA_low_plot <- ggplot(PVA_low_all_long) +
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
  ylim(0, 100) +
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
PVA_low_plot




PVA_middle_4clutch <- as.data.frame(cbind(
  years,
  PVA_min,
  PVA_med,
  PVA_max
))

PVA_middle_6clutch <- as.data.frame(cbind(
  years,
  PVA_min,
  PVA_med,
  PVA_max
))

PVA_middle_10clutch <- as.data.frame(cbind(
  years,
  PVA_min,
  PVA_med,
  PVA_max
))

PVA_middle_allclutch <- rbind(PVA_middle_4clutch,
                             PVA_middle_6clutch,
                             PVA_middle_10clutch)
PVA_middle_allclutch$clutchHS <- rep(c(4, 6, 10), each = 31)

names(PVA_middle_allclutch) <- c("years",
                                "low values",
                                "med values",
                                "high values",
                                "clutches HS")

PVA_upper_4clutch <- as.data.frame(cbind(
  years,
  PVA_min,
  PVA_med,
  PVA_max
))

PVA_upper_6clutch <- as.data.frame(cbind(
  years,
  PVA_min,
  PVA_med,
  PVA_max
))

PVA_upper_10clutch <- as.data.frame(cbind(
  years,
  PVA_min,
  PVA_med,
  PVA_max
))

PVA_upper_allclutch <- rbind(PVA_upper_4clutch,
                             PVA_upper_6clutch,
                             PVA_upper_10clutch)
PVA_upper_allclutch$clutchHS <- rep(c(4, 6, 10), each = 31)

names(PVA_upper_allclutch) <- c("years",
                                "low values",
                                "med values",
                                "high values",
                                "clutches HS")


PVA_20_4clutch <- as.data.frame(cbind(
  years,
  PVA_min,
  PVA_med,
  PVA_max
))

PVA_20_6clutch <- as.data.frame(cbind(
  years,
  PVA_min,
  PVA_med,
  PVA_max
))

PVA_20_10clutch <- as.data.frame(cbind(
  years,
  PVA_min,
  PVA_med,
  PVA_max
))

PVA_lower_allclutch <- rbind(PVA_lower_4clutch,
                             PVA_lower_6clutch,
                             PVA_lower_10clutch)
PVA_lower_allclutch$clutchHS <- rep(c(4, 6, 10), each = 31)

names(PVA_lower_allclutch) <- c("years",
                                "low values",
                                "med values",
                                "high values",
                                "clutches HS")


names(PVA_glen) <- c("years",
                     "low values",
                     "med values",
                     "high values")

names(PVA_lower) <- c("years",
                      "low values",
                      "med values",
                      "high values")

names(PVA_middle) <- c("years",
                       "low values",
                       "med values",
                       "high values")

names(PVA_upper) <- c("years",
                      "low values",
                      "med values",
                      "high values")

names(PVA_20) <- c("years",
                   "low values",
                   "med values",
                   "high values")



siteRemPVA <- as_tibble(rbind(
  PVA_lower,
  PVA_middle,
  PVA_upper,
  PVA_20
))


names(PVAoutput) <- c("years",
                      "low values",
                      "med values",
                      "high values")

values <- "values"
pop_size <- "population size (individs)"
gathercols <- c("low values", "med values", "high values")

PVAoutput.long <- gather(PVAoutput, values, pop_size, gathercols)

HS_10_long <- gather(HS_10_output, values, pop_size, gathercols)
HS_15_long <- gather(HS_15_output, values, pop_size, gathercols)
HS_20_long <- gather(HS_20_output, values, pop_size, gathercols)
HS_25_long <- gather(HS_25_output, values, pop_size, gathercols)
HS_30_long <- gather(HS_30_output, values, pop_size, gathercols)

HS_35_long <- gather(HS_35_output, values, pop_size, gathercols)
a <- rbind(HS_10_long,
           HS_15_long,
           HS_20_long,
           HS_25_long,
           HS_30_long,
           HS_35_long)
a$clutches <- rep(c(10, 15, 20, 25, 30, 35), each  = 93)

names(HS_10_output) <- c("years",
                         "low values",
                         "med values",
                         "high values")

names(HS_100_output) <- c("years",
                          "low values",
                          "med values",
                          "high values")

# Plot
PVAoutput.plot <- ggplot(PVAoutput.long) +
  geom_line(aes(x = years,
                y = pop_size,
                group = values,
                colour = values)) +
  xlab("Years") +
  ylab("Curlew population size (individuals)") +
  ylim(0, 800) +
  scale_x_continuous(breaks = seq(2022, 2052, 1)) +
  theme_bw() + 
  theme(axis.title = element_text(size = 32),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        legend.position = "none")
PVAoutput.plot


PVAoutput.plot <- ggplot(a) +
  facet_wrap(~clutches) + 
  geom_line(aes(x = years,
                y = pop_size,
                group = values,
                colour = values)) +
  stat_summary(aes(y = pop_size,x = years, group = clutches),
               geom = "smooth",
               fun = "mean",
               colour = "black") +
  xlab("Years") +
  ylab("Curlew population size (individuals)") +
  ylim(0, 1000) +
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
PVAoutput.plot


NI.sim.out.plot.ribbon <- NI.sim.out.plot +
  geom_ribbon(aes(x = years,
                  ymin = lcl,
                  ymax = ucl),
              alpha = 0.3) +
  xlab("Years") +
  ylab("Curlew population size (individuals)") +
  ylim(0, 1500) +
  theme_bw() + 
  theme(axis.title = element_text(size = 32),
        axis.text = element_text(size = 16, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        legend.position = "none")










#########################################################


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
        "data/NI_HS_A3_sims_hphr_1.RDS")



#----- Population plots -----
NI.sim.out <- readRDS("data/NI_HS_base_A1_sims_hplr_1.RDS")
names.NI.sim.out <- rep(as.character(c(NI.pop.size)), (56000/5))

names(NI.sim.out) <- names.NI.sim.out
NI.sim.out <- sim.out

years <- seq(2022, 2052, 1)

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


NI.sim.out.df.long <- data.frame(years = rep(years, times = 5),
                                 initN = rep(NI.pop.size, each = 31),
                                 mean = c(NI.sim.out.200.mean,
                                          NI.sim.out.225.mean,
                                          NI.sim.out.250.mean,
                                          NI.sim.out.275.mean,
                                          NI.sim.out.300.mean),
                                 lcl = c(NI.sim.out.200.lcl,
                                         NI.sim.out.225.lcl,
                                         NI.sim.out.250.lcl,
                                         NI.sim.out.275.lcl,
                                         NI.sim.out.300.lcl),
                                 ucl = c(NI.sim.out.200.ucl,
                                         NI.sim.out.225.ucl,
                                         NI.sim.out.250.ucl,
                                         NI.sim.out.275.ucl,
                                         NI.sim.out.300.ucl)
)


write.table(NI.sim.out.df.long,
            "data/NI_HS_A3_means_hphr_fullsims_1.csv",
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
  ylim(0, 1500) +
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
bettermatch.200 <- which(NI.sim.out.200[31 ,] < 10,
                         arr.ind = TRUE)
NI.sim.out.200.match <- NI.sim.out.200[, bettermatch.200]
NI.sim.out.200.match.mean <- rowMeans(NI.sim.out.200.match)
NI.sim.out.200.match.lcl <- apply(NI.sim.out.200.match, 1,
                                  quantile, probs = 0.1)
NI.sim.out.200.match.ucl <- apply(NI.sim.out.200.match, 1,
                                  quantile, probs = 0.9)

bettermatch.225 <- which(NI.sim.out.225[31 ,] < 10,
                         arr.ind = TRUE)
NI.sim.out.225.match <- NI.sim.out.225[, bettermatch.225]
NI.sim.out.225.match.mean <- rowMeans(NI.sim.out.225.match)
NI.sim.out.225.match.lcl <- apply(NI.sim.out.225.match, 1,
                                  quantile, probs = 0.1)
NI.sim.out.225.match.ucl <- apply(NI.sim.out.225.match, 1,
                                  quantile, probs = 0.9)

bettermatch.250 <- which(NI.sim.out.250[31 ,] < 10,
                         arr.ind = TRUE)
NI.sim.out.250.match <- NI.sim.out.250[, bettermatch.250]
NI.sim.out.250.match.mean <- rowMeans(NI.sim.out.250.match)
NI.sim.out.250.match.lcl <- apply(NI.sim.out.250.match, 1,
                                  quantile, probs = 0.1)
NI.sim.out.250.match.ucl <- apply(NI.sim.out.250.match, 1,
                                  quantile, probs = 0.9)

bettermatch.275 <- which(NI.sim.out.275[31 ,] < 10,
                         arr.ind = TRUE)
NI.sim.out.275.match <- NI.sim.out.275[, bettermatch.275]
NI.sim.out.275.match.mean <- rowMeans(NI.sim.out.275.match)
NI.sim.out.275.match.lcl <- apply(NI.sim.out.275.match, 1,
                                  quantile, probs = 0.1)
NI.sim.out.275.match.ucl <- apply(NI.sim.out.275.match, 1,
                                  quantile, probs = 0.9)

bettermatch.300 <- which(NI.sim.out.300[31 ,] < 10,
                         arr.ind = TRUE)
NI.sim.out.300.match <- NI.sim.out.300[, bettermatch.300]
NI.sim.out.300.match.mean <- rowMeans(NI.sim.out.300.match)
NI.sim.out.300.match.lcl <- apply(NI.sim.out.300.match, 1,
                                  quantile, probs = 0.1)
NI.sim.out.300.match.ucl <- apply(NI.sim.out.300.match, 1,
                                  quantile, probs = 0.9)


NI.sim.out.df.match <- data.frame(years = rep(years, times = 5),
                                  initN = rep(NI.pop.size, each = 31),
                                  mean = c(NI.sim.out.200.match.mean,
                                           NI.sim.out.225.match.mean,
                                           NI.sim.out.250.match.mean,
                                           NI.sim.out.275.match.mean,
                                           NI.sim.out.300.match.mean),
                                  lcl = c(NI.sim.out.200.match.lcl,
                                          NI.sim.out.225.match.lcl,
                                          NI.sim.out.250.match.lcl,
                                          NI.sim.out.275.match.lcl,
                                          NI.sim.out.300.match.lcl),
                                  ucl = c(NI.sim.out.200.match.ucl,
                                          NI.sim.out.225.match.ucl,
                                          NI.sim.out.250.match.ucl,
                                          NI.sim.out.275.match.ucl,
                                          NI.sim.out.300.match.ucl)
)


write.table(NI.sim.out.df.match,
            "data/NI_HS_A3_hphr_means_match10_1.csv",
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
  ylim(0, 1000) +
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



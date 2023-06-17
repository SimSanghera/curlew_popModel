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

# Update 03/02/22:
#   using the same stages (excluding headstarted birds) as the base and 
#   scenario model to compare to the two stage. 
#   This model will also be split into low philopatry and high philopatry 
#   cases, a) to reduce processing time, b) compare between these params.


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
# Update: 03/02/2022 - using reduced number of values for parameters
#   with all stages as would be used in the base and scenario models.
#   for comparison with previously run historic model.
#   no HS birds, so no need for intro.years

# Pop Sizes NI
pop.1987 <- seq(3800, 6250, 400)
pop.1999 <- seq(1240, 2928, 200)
pop.2013 <- seq(252, 783, 200)

# Pop size ROI
ROI.pop.1985 <- seq(5000, 9000, 1000)


#-----  Parameters  -----

# Population Size
#NI.pop.size   <- seq(200, 300, 25)      # breeding pairs

#ROI.pop.size  <- seq(100, 200, 25)


# Survival Rates & Curlew Related Parameters
S.ad          <- seq(0.78, 0.92, 0.04)  # Cook et al., 2021
# Roodbergen et al., 2021

S.jw          <- seq(0.20, 0.40, 0.05)  # Wild bird survival
# Cook et al., 2021

S.jhs         <- seq(0.20, 0.40, 0.05)  # HS bird survival
# Assume same as wild for now
# Geoff Hilton (pers comms)

S.sub         <- seq(0.78, 0.92, 0.04)  # subadult survival

S.imm         <- seq(0.40, 0.70, 0.1)  # survival of immature (age fledge - 1)

prod.rates    <- seq(0.14, 0.52, 0.06)  # Grant et al., 1999
# Harris et al., 2019
# Zielonka et al., 2019
# Irish Breeding Bird Report

hs.pop        <- 0                      # baseline - 0 HS

capt.fail     <- 0.2                    # The number of HS birds that
# fail during rearing & can't be
# released

intro.years   <- 1           # years in which hs birds will be
# introduced

phil.p        <- 0.3                    # probability of first years return
# to release site

return.p      <- 0.90                   # probability of breeding birds 
# returning to release site


# Environmental Parameters
NI.K          <- 15000                  # carrying capacity of NI

ROI.K         <- 20000                  # carrying capacity of ROI

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


#-----  Create Combinations of Demographic Parameters -----

simul.in  <- expand.grid(init.N       = pop.1987,
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
            "data/NI_historic_4stage_hphr_parameters_2.csv",
            sep = ",",
            row.names = FALSE)



#-----  Defining Functions for Stochastic PVA -----

# Defining the structure of the population model with adult & juvenile stages
# Need to add other age stages representing survival
#   - egg to hatch
#   - 1 year breeders
#   - S.jw (and S.jr) is the same as prod.rate as we are assuming same values
# And prod.rate is the proportion of chicks that make it to next stage

bird.matrix <- expression(
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
        pop.size  <- c(0, hs.numbers, 0, 0, (init.N * 2))
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
      if (y < hs.years) next.year <- next.year + hs.numbers
      
      # Put in array
      Pop.Array[y, rep] <- next.year
      
    } # end year
    
  }   # end rep
  
  return(Pop.Array)
  
}



#----- Test -----
s = 1
rep = 1
y = 2
init.N           = simul.in[s, 1]
F.loop           = simul.in[s, 7]
S.ad.loop        = simul.in[s, 2]
S.jw.loop        = simul.in[s, 3]
S.jhs.loop       = simul.in[s, 4]
S.imm.loop       = simul.in[s, 5]
S.sub.loop       = simul.in[s, 6]
philopatry       = simul.in[s, 8]
sitefidelity     = simul.in[s, 9]
hs.numbers       = simul.in[s, 10]
hs.years         = simul.in[s, 11]


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
####################################



#-----  Run Sims with Parallel Processing -----

# Boost speed using parallel processing
cl <- makeCluster(8)
registerDoParallel(cl, cores = 8)

# Start loop
system.time(
  
  sim.out <- foreach(s = c(1:dim(simul.in)[1]),
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






names.NI.sim.out <- rep(as.character(c(pop.1987)), (78400/7))




# Write to save file
saveRDS(sim.out,
        "data/NI_historic_4stage_hphr_sims_1.RDS")
# lphr = phil 0.3, ret = 0.9, breed.fail = 0.2
# hphr = phil 0.7, ret = 0.9, breed.fail = 0.2


#-----  Plotting  -----
NI.sim.out <- readRDS("data/NI_historic_4stage_hphr_sims_1.RDS")

names(NI.sim.out) <- names.NI.sim.out

years <- seq(1987, 2017, 1)

NI.sim.out.3800 <- NI.sim.out[names(NI.sim.out) == "3800"]
NI.sim.out.3800 <- do.call(cbind, NI.sim.out.3800)
NI.sim.out.3800.mean <- rowMeans(NI.sim.out.3800)
NI.sim.out.3800.lcl  <- apply(NI.sim.out.3800, 1, quantile, probs = 0.1)
NI.sim.out.3800.ucl  <- apply(NI.sim.out.3800, 1, quantile, probs = 0.9)

NI.sim.out.4200 <- NI.sim.out[names(NI.sim.out) == "4200"]
NI.sim.out.4200 <- do.call(cbind, NI.sim.out.4200)
NI.sim.out.4200.mean <- rowMeans(NI.sim.out.4200)
NI.sim.out.4200.lcl  <- apply(NI.sim.out.4200, 1, quantile, probs = 0.1)
NI.sim.out.4200.ucl  <- apply(NI.sim.out.4200, 1, quantile, probs = 0.9)

NI.sim.out.4600 <- NI.sim.out[names(NI.sim.out) == "4600"]
NI.sim.out.4600 <- do.call(cbind, NI.sim.out.4600)
NI.sim.out.4600.mean <- rowMeans(NI.sim.out.4600)
NI.sim.out.4600.lcl  <- apply(NI.sim.out.4600, 1, quantile, probs = 0.1)
NI.sim.out.4600.ucl  <- apply(NI.sim.out.4600, 1, quantile, probs = 0.9)

NI.sim.out.5000 <- NI.sim.out[names(NI.sim.out) == "5000"]
NI.sim.out.5000 <- do.call(cbind, NI.sim.out.5000)
NI.sim.out.5000.mean <- rowMeans(NI.sim.out.5000)
NI.sim.out.5000.lcl  <- apply(NI.sim.out.5000, 1, quantile, probs = 0.1)
NI.sim.out.5000.ucl  <- apply(NI.sim.out.5000, 1, quantile, probs = 0.9)

NI.sim.out.5400 <- NI.sim.out[names(NI.sim.out) == "5400"]
NI.sim.out.5400 <- do.call(cbind, NI.sim.out.5400)
NI.sim.out.5400.mean <- rowMeans(NI.sim.out.5400)
NI.sim.out.5400.lcl  <- apply(NI.sim.out.5400, 1, quantile, probs = 0.1)
NI.sim.out.5400.ucl  <- apply(NI.sim.out.5400, 1, quantile, probs = 0.9)

NI.sim.out.5800 <- NI.sim.out[names(NI.sim.out) == "5800"]
NI.sim.out.5800 <- do.call(cbind, NI.sim.out.5800)
NI.sim.out.5800.mean <- rowMeans(NI.sim.out.5800)
NI.sim.out.5800.lcl  <- apply(NI.sim.out.5800, 1, quantile, probs = 0.1)
NI.sim.out.5800.ucl  <- apply(NI.sim.out.5800, 1, quantile, probs = 0.9)

NI.sim.out.6200 <- NI.sim.out[names(NI.sim.out) == "6200"]
NI.sim.out.6200 <- do.call(cbind, NI.sim.out.6200)
NI.sim.out.6200.mean <- rowMeans(NI.sim.out.6200)
NI.sim.out.6200.lcl  <- apply(NI.sim.out.6200, 1, quantile, probs = 0.1)
NI.sim.out.6200.ucl  <- apply(NI.sim.out.6200, 1, quantile, probs = 0.9)


NI.sim.out.df.long <- data.frame(years = rep(years, times = 7),
                             initN = rep(pop.1987, each = 31),
                             mean = c(NI.sim.out.3800.mean,
                                      NI.sim.out.4200.mean,
                                      NI.sim.out.4600.mean,
                                      NI.sim.out.5000.mean,
                                      NI.sim.out.5400.mean,
                                      NI.sim.out.5800.mean,
                                      NI.sim.out.6200.mean),
                             lcl = c(NI.sim.out.3800.lcl,
                                     NI.sim.out.4200.lcl,
                                     NI.sim.out.4600.lcl,
                                     NI.sim.out.5000.lcl,
                                     NI.sim.out.5400.lcl,
                                     NI.sim.out.5800.lcl,
                                     NI.sim.out.6200.lcl),
                             ucl = c(NI.sim.out.3800.ucl,
                                     NI.sim.out.4200.ucl,
                                     NI.sim.out.4600.ucl,
                                     NI.sim.out.5000.ucl,
                                     NI.sim.out.5400.ucl,
                                     NI.sim.out.5800.ucl,
                                     NI.sim.out.6200.ucl)
                             )


write.table(NI.sim.out.df.long,
            "data/NI_sim_out_historic_plot_1987_hphr_1.csv",
            sep = ",",
            row.names = FALSE)



# Create points for 3 survey years
pop <- c((2*5000), (2*2091), (2*526))
pop.lcl <- c((2*3800), (2*1243), (2*252))
pop.ucl <- c((2*6250), (2*2928), (2*783))
NI.pop.ests <- data.frame(years = c(1987, 1999, 2013),
                          popsize = pop,
                          pop.min = pop.lcl,
                          pop.max = pop.ucl)


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
                 y = popsize),
             colour = "black") +
  geom_errorbar(data = NI.pop.ests,
                aes(x = years, ymin = pop.min, ymax = pop.max),
                colour = "black",
                width = 1,
                size = 0.5) +
  xlab("Years") + 
  ylab("Curlew Population Size (individuals)") +
  ylim(0, 20000) + 
  theme(panel.background = element_rect(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = "none")

NI.sim.out.plot
NI.sim.out.plot.ribbon

# Extra theme edits
# panel.background = element_rect(),
# axis.title = element_text(size = 32),
# axis.text = element_text(size = 16, color = "black"),
# panel.grid.major = element_blank(),
# panel.grid.minor = element_blank(),
# panel.border = element_blank(),
# legend.position = "none",
# legend.background = element_rect(fill = "transparent"),
# legend.box.background = element_rect(fill = "transparent"),
# legend.text = element_text(size = 28),
# legend.title = element_text(size = 28)


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


bettermatch <- which(NI.sim.out.4200[31 ,] < 500,
                     arr.ind = TRUE)
NI.sim.out.4200.match <- NI.sim.out.4200[, bettermatch]
NI.sim.out.4200.match.mean <- rowMeans(NI.sim.out.4200.match)
NI.sim.out.4200.match.lcl <- apply(NI.sim.out.4200.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.4200.match.ucl <- apply(NI.sim.out.4200.match, 1,
                                   quantile, probs = 0.9)


bettermatch <- which(NI.sim.out.4600[31 ,] < 500,
                     arr.ind = TRUE)
NI.sim.out.4600.match <- NI.sim.out.4600[, bettermatch]
NI.sim.out.4600.match.mean <- rowMeans(NI.sim.out.4600.match)
NI.sim.out.4600.match.lcl <- apply(NI.sim.out.4600.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.4600.match.ucl <- apply(NI.sim.out.4600.match, 1,
                                   quantile, probs = 0.9)


bettermatch <- which(NI.sim.out.5000[31 ,] < 500,
                     arr.ind = TRUE)
NI.sim.out.5000.match <- NI.sim.out.5000[, bettermatch]
NI.sim.out.5000.match.mean <- rowMeans(NI.sim.out.5000.match)
NI.sim.out.5000.match.lcl <- apply(NI.sim.out.5000.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.5000.match.ucl <- apply(NI.sim.out.5000.match, 1,
                                   quantile, probs = 0.9)


bettermatch <- which(NI.sim.out.5400[31 ,] < 500,
                     arr.ind = TRUE)
NI.sim.out.5400.match <- NI.sim.out.5400[, bettermatch]
NI.sim.out.5400.match.mean <- rowMeans(NI.sim.out.5400.match)
NI.sim.out.5400.match.lcl <- apply(NI.sim.out.5400.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.5400.match.ucl <- apply(NI.sim.out.5400.match, 1,
                                   quantile, probs = 0.9)


bettermatch <- which(NI.sim.out.5800[31 ,] < 500,
                     arr.ind = TRUE)
NI.sim.out.5800.match <- NI.sim.out.5800[, bettermatch]
NI.sim.out.5800.match.mean <- rowMeans(NI.sim.out.5800.match)
NI.sim.out.5800.match.lcl <- apply(NI.sim.out.5800.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.5800.match.ucl <- apply(NI.sim.out.5800.match, 1,
                                   quantile, probs = 0.9)


bettermatch <- which(NI.sim.out.6200[31 ,] < 500,
                     arr.ind = TRUE)
NI.sim.out.6200.match <- NI.sim.out.6200[, bettermatch]
NI.sim.out.6200.match.mean <- rowMeans(NI.sim.out.6200.match)
NI.sim.out.6200.match.lcl <- apply(NI.sim.out.6200.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.6200.match.ucl <- apply(NI.sim.out.6200.match, 1,
                                   quantile, probs = 0.9)


NI.sim.out.df.match <- data.frame(years = rep(years, times = 7),
                                 initN = rep(pop.1987, each = 31),
                                 mean = c(NI.sim.out.3800.match.mean,
                                          NI.sim.out.4200.match.mean,
                                          NI.sim.out.4600.match.mean,
                                          NI.sim.out.5000.match.mean,
                                          NI.sim.out.5400.match.mean,
                                          NI.sim.out.5800.match.mean,
                                          NI.sim.out.6200.match.mean),
                                 lcl = c(NI.sim.out.3800.match.lcl,
                                         NI.sim.out.4200.match.lcl,
                                         NI.sim.out.4600.match.lcl,
                                         NI.sim.out.5000.match.lcl,
                                         NI.sim.out.5400.match.lcl,
                                         NI.sim.out.5800.match.lcl,
                                         NI.sim.out.6200.match.lcl),
                                 ucl = c(NI.sim.out.3800.match.ucl,
                                         NI.sim.out.4200.match.ucl,
                                         NI.sim.out.4600.match.ucl,
                                         NI.sim.out.5000.match.ucl,
                                         NI.sim.out.5400.match.ucl,
                                         NI.sim.out.5800.match.ucl,
                                         NI.sim.out.6200.match.ucl)
)

write.table(NI.sim.out.df.match,
            "data/NI_sim_out_historic_plot_1987_hphr_match_1.csv",
            sep = ",",
            row.names = FALSE)

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
                 y = popsize),
             colour = "black") +
  geom_errorbar(data = NI.pop.ests,
                aes(x = years, ymin = pop.min, ymax = pop.max),
                colour = "black",
                width = 1,
                size = 0.5) +
  xlab("Years") + 
  ylab("Curlew Population Size (individuals)") +
  ylim(0, 20000) + 
  theme(panel.background = element_rect(),
        axis.title = element_text(size = 32),
        axis.text = element_text(size = 16, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = "none")

NI.pop.ests <- data.frame(years = c(1987, 1999, 2013),
                          popsize = c((2*5000), (2*2091), (2*526)))



# Get combinations of parameters
bettermatch.3800 <- which(NI.sim.out.3800[31 ,] < 500,
                     arr.ind = TRUE)

NI.3800.match.params <- NI.sim.params[bettermatch.3800 ,]


##############################################################
##############################################################
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
        pop.size  <- c(0, hs.numbers, 0, 0, (init.N * 2))
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
                         trunc(ROI.Ricker(Pop.Array[y - 1, rep])))
      
      # Catastrophe
      if(runif(1) < catastrophe.p) next.year <- next.year * catastrophe.i
      if (y < hs.years) next.year <- next.year + hs.numbers
      
      # Put in array
      Pop.Array[y, rep] <- next.year
      
    } # end year
    
  }   # end rep
  
  return(Pop.Array)
  
}



# Combination of parameters
ROI.simul.in <- expand.grid(init.N       = ROI.pop.1985,
                            S.ad.loop    = S.ad,
                            S.jw.loop    = S.jw,
                            S.jhs.loop   = S.jhs,
                            S.imm.loop   = S.imm,
                            S.sub.loop   = S.sub,
                            F.loop       = prod.rates,
                            philopatry   = phil.p,
                            sitefidelity = return.p,
                            hs.numbers   = hs.pop,
                            hs.years     = intro.years) # for release, add intro years

dim(simul.in)
names(ROI.simul.in) <- c("initial.pop",
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


write.table(ROI.simul.in, "data/ROI_historic_4stage_hphr_params_1.csv",
            sep = ",",
            row.names = FALSE)


#-----  Try Parallel Processing -----


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
                       init.N <- ROI.simul.in[s,1]
                       default <- ROI.PVA(
                         n.reps,
                         n.years,
                         init.N,
                         ROI.K,
                         catastrophe.p,
                         catastrophe.i,
                         F.loop           = ROI.simul.in[s, 7],
                         S.ad.loop        = ROI.simul.in[s, 2],
                         S.jw.loop        = ROI.simul.in[s, 3],
                         S.jhs.loop       = ROI.simul.in[s, 4],
                         S.imm.loop       = ROI.simul.in[s, 5],
                         S.sub.loop       = ROI.simul.in[s, 6],
                         philopatry       = ROI.simul.in[s, 8],
                         sitefidelity     = ROI.simul.in[s, 9],
                         hs.numbers       = ROI.simul.in[s, 10],
                         hs.years         = ROI.simul.in[s, 11],
                         capt.fail,
                         breed.fail
                       )
                       
                     }
  
)

# Stop parallel processing
stopCluster(cl)

roi.names <- rep(as.character(ROI.pop.1985), (56000/5))
names(ROI.sim.out) <- roi.names


saveRDS(ROI.sim.out,
        "data/ROI_historic_4stage_lphr_sims_1.RDS")
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



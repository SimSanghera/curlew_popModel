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
pop.1987 <- seq(3800, 6250, 550)
pop.1999 <- c(1243, 2091, 2928)
pop.2013 <- c(252, 526, 783)

# Pop size ROI
ROI.pop.1985 <- seq(6000, 8000, 500)


S.ad            <- c(0.78, 0.86, 0.92)  # adult survival
# rates obtained from literature
# using UK wide rates
# NEED TO OBTAIN EU RATES

S.jw            <- seq(0.2, 0.40, 0.05)  # Juvenile survival - chick to
# fledging. UK wide. GET EU Numbers

S.jhs         <- seq(0.20, 0.40, 0.05)  # HS bird survival
# Assume same as wild for now
# Geoff Hilton (pers comms)

S.sub         <- c(0.78, 0.86, 0.92)  # subadult survival

S.imm         <- seq(0.40, 0.70, 0.15)  # survival of immature (age fledge - 1)

prod.rates    <- seq(0.14, 0.52, 0.06)  # Grant et al., 1999
# Harris et al., 2019
# Zielonka et al., 2019
# Irish Breeding Bird Report

hs.clutches   <- 0
hs.pop        <- (4* hs.clutches)       # baseline - 0 HS

capt.fail     <- 0.2                    # The number of HS birds that
# fail during rearing & can't be
# released

intro.years   <- 5           # years in which hs birds will be
# introduced

phil.p        <- 0.3                    # probability of first years return
# to release site

return.p      <- 0.70                   # probability of breeding birds 
# returning to release site


# Environmental Parameters
NI.K          <- 25000                  # carrying capacity of NI

ROI.K         <- 25000                  # carrying capacity of ROI

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
  # f_w           imm                     sub                     adult  
  0,              (prod.rates*0.5)*0.005, (prod.rates*0.5)*0.1, prod.rates*0.5,
  S.jw*phil.p,      0,                      0,                      0,
  0,                S.imm*return.p,         0,                      0,
  0,                0,                      S.sub*return.p,     S.ad*return.p
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



#-----  Create Combinations of Demographic Parameters -----

simul.in  <- expand.grid(init.N       = ROI.pop.1985, # set as starting pop
                         S.ad.loop    = S.ad,
                         S.jw.loop    = S.jw,
                         S.imm.loop   = S.imm,
                         S.sub.loop   = S.sub,
                         F.loop       = prod.rates,
                         philopatry   = phil.p,
                         sitefidelity = return.p)

names(simul.in) <- c("initial.pop",
                     "S.ad",
                     "S.jw",
                     "S.imm",
                     "S.sub",
                     "Prod.Rates",
                     "philopatry",
                     "site.fidelity")

write.table(simul.in,
            "data/curlew_NI_historic_parameters_2.csv",
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
  S.imm.loop,
  S.sub.loop,
  philopatry,
  sitefidelity,
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
      bird.vr <- list(
                      init.N = init.N,
                      prod.rates  = F.Year,
                      S.ad        = S.ad.loop,
                      S.jw        = S.jw.loop,
                      S.imm       = S.imm.loop,
                      S.sub       = S.sub.loop,
                      phil.p      = philopatry,
                      return.p    = sitefidelity
      )
      
      # Create leslie matrix
      a <- matrix(sapply(bird.matrix,
                         eval,
                         bird.vr,
                         NULL),
                  nrow = sqrt(length(bird.matrix)),
                  byrow = TRUE)
      
      # Starting pop vector
      # First year has no wild chicks, as pop counts only report adults
      # Year 2 will be any chicks born to the initial adults
      #                     0 imm, 0 sub
      #                     surviving & returning adults
      # year 3 onwards will be stages from pop.projection

      
      if (y == 2) {
        pop.size  <- c(((F.loop * 0.5) * (init.N*2)), # prod per pair * individs
                       0, # no imm
                       0, # no sub
                       (init.N*2)*S.ad.loop*sitefidelity) # surv & return adults

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
# 
# init.N          = simul.in[s, 1]
# F.loop          = simul.in[s, 6]
# S.ad.loop       = simul.in[s, 2]
# S.jw.loop       = simul.in[s, 3]
# S.imm.loop      = simul.in[s, 4]
# S.sub.loop      = simul.in[s, 5]
# philopatry      = simul.in[s, 7]
# sitefidelity    = simul.in[s, 8]


#----- Run Sims with Parallel Processing  -----
cl <- makeCluster(8)
registerDoParallel(cl, cores = 8)

# start loop
system.time(
  
  ROI.sim.out <- foreach(s = c(1:dim(simul.in)[1]),
                        .packages = "popbio") %dopar% {
                          
                          # start PVA
                          init.N            = simul.in[s, 1]
                          default <- PVA(
                            n.reps,
                            n.years,
                            init.N,
                            K,
                            catastrophe.p,
                            catastrophe.i,
                            F.loop          = simul.in[s, 6],
                            S.ad.loop       = simul.in[s, 2],
                            S.jw.loop       = simul.in[s, 3],
                            S.imm.loop      = simul.in[s, 4],
                            S.sub.loop      = simul.in[s, 5],
                            philopatry       = simul.in[s, 7],
                            sitefidelity    = simul.in[s, 8],
                            breed.fail
                          )
                          
                        }
  
)

# Stop processing
stopCluster(cl)



testPVA <- PVA(
  init.N          = simul.in[s, 1],
  n.reps,
  n.years,
  init.N,
  K,
  catastrophe.p,
  catastrophe.i,
  F.loop          = simul.in[s, 6],
  S.ad.loop       = simul.in[s, 2],
  S.jw.loop       = simul.in[s, 3],
  S.imm.loop      = simul.in[s, 4],
  S.sub.loop      = simul.in[s, 5],
  philopatry      = simul.in[s, 7],
  sitefidelity    = simul.in[s, 8],
  breed.fail
)

testPVA.means <- as.data.frame(rowMeans(testPVA))
testPVA.lcl   <- apply(testPVA, 1, quantile, probs = 0.1)
testPVA.ucl   <- apply(testPVA, 1, quantile, probs = 0.9)



#----- Test ------
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






# Write to save file
saveRDS(NI.sim.out,
        "data/NI_historic_mod_1987_sims_5.RDS")
#ver 3 = breed.fail 0,4

#-----  Plotting  -----
NI.sim.out.hist <- readRDS("data/NI_historic_mod_1987_sims_5.RDS")
ROI.sim.out.hist <- readRDS("data/ROI_historic1985_sims.RDS")

names.NI.sim.out.hist <- rep(as.character(c(pop.1987)), (4725/5))

names(NI.sim.out.hist) <- names.NI.sim.out.hist

years.hist <- seq(1987, 2017, 1)

# Extract sims for each starting pop
NI.sim.out.3800 <- NI.sim.out.hist[names(NI.sim.out.hist) == "3800"]
NI.sim.out.3800 <- do.call(cbind, NI.sim.out.3800)
NI.sim.out.3800.mean <- rowMeans(NI.sim.out.3800)
NI.sim.out.3800.lcl  <- apply(NI.sim.out.3800, 1, quantile, probs = 0.1)
NI.sim.out.3800.ucl  <- apply(NI.sim.out.3800, 1, quantile, probs = 0.9)

NI.sim.out.4350 <- NI.sim.out.hist[names(NI.sim.out.hist) == "4350"]
NI.sim.out.4350 <- do.call(cbind, NI.sim.out.4350)
NI.sim.out.4350.mean <- rowMeans(NI.sim.out.4350)
NI.sim.out.4350.lcl  <- apply(NI.sim.out.4350, 1, quantile, probs = 0.1)
NI.sim.out.4350.ucl  <- apply(NI.sim.out.4350, 1, quantile, probs = 0.9)

NI.sim.out.4900 <- NI.sim.out.hist[names(NI.sim.out.hist) == "4900"]
NI.sim.out.4900 <- do.call(cbind, NI.sim.out.4900)
NI.sim.out.4900.mean <- rowMeans(NI.sim.out.4900)
NI.sim.out.4900.lcl  <- apply(NI.sim.out.4900, 1, quantile, probs = 0.1)
NI.sim.out.4900.ucl  <- apply(NI.sim.out.4900, 1, quantile, probs = 0.9)

NI.sim.out.5450 <- NI.sim.out.hist[names(NI.sim.out.hist) == "5450"]
NI.sim.out.5450 <- do.call(cbind, NI.sim.out.5450)
NI.sim.out.5450.mean <- rowMeans(NI.sim.out.5450)
NI.sim.out.5450.lcl  <- apply(NI.sim.out.5450, 1, quantile, probs = 0.1)
NI.sim.out.5450.ucl  <- apply(NI.sim.out.5450, 1, quantile, probs = 0.9)

# remove nan
# NI.sim.out.6000 <- NI.sim.out.6000[, colSums(is.na(NI.sim.out.6000)) == 0]

NI.sim.out.6000 <- NI.sim.out.hist[names(NI.sim.out.hist) == "6000"]
NI.sim.out.6000 <- do.call(cbind, NI.sim.out.6000)
NI.sim.out.6000.mean <- rowMeans(NI.sim.out.6000)
NI.sim.out.6000.lcl  <- apply(NI.sim.out.6000, 1, quantile, probs = 0.1)
NI.sim.out.6000.ucl  <- apply(NI.sim.out.6000, 1, quantile, probs = 0.9)



# Combine into tibble

NI.sim.out.df <- tibble(years = as.factor(years.hist), 
                         NI.sim.out.3800.mean,
                         NI.sim.out.3800.lcl,
                         NI.sim.out.3800.ucl,
                         NI.sim.out.4350.mean,
                         NI.sim.out.4350.lcl,
                         NI.sim.out.4350.ucl,
                         NI.sim.out.4900.mean,
                         NI.sim.out.4900.lcl,
                         NI.sim.out.4900.ucl,
                         NI.sim.out.5450.mean,
                         NI.sim.out.5450.lcl,
                         NI.sim.out.5450.ucl,
                         NI.sim.out.6000.mean,
                         NI.sim.out.6000.lcl,
                         NI.sim.out.6000.ucl)

write.table(I.sim.out.df,
            "data/NI_historic_1987_df_1.csv",
            sep = ",",
            row.names = FALSE)

NI.init.N <- rep(pop.1987, 3)

NI.sim.out.df.long <- data.frame(years = rep(years, times = 5),
                             initN = rep(pop.1987, each = 31),
                             mean = c(NI.sim.out.3800.mean,
                                          NI.sim.out.4350.mean,
                                          NI.sim.out.4900.mean,
                                          NI.sim.out.5450.mean,
                                          NI.sim.out.6000.mean),
                             lcl = c(NI.sim.out.3800.lcl,
                                         NI.sim.out.4350.lcl,
                                         NI.sim.out.4900.lcl,
                                         NI.sim.out.5450.lcl,
                                         NI.sim.out.6000.lcl),
                             ucl = c(NI.sim.out.3800.ucl,
                                         NI.sim.out.4350.ucl,
                                         NI.sim.out.4900.ucl,
                                         NI.sim.out.5450.ucl,
                                         NI.sim.out.6000.ucl)
                             )


write.table(NI.sim.out.df.long,
            "data/NI_sim_out_historic_1987.csv",
            sep = ",",
            row.names = FALSE)

NIpops <- tibble(years = c("1987", "1999", "2013"),
                 popsize = c(5000, 2091, 526),
                 popsize_l = c(3800, 1243, 252),
                 popsize_u = c(6250, 2928, 783))



kispal <- fish(n = 8, option = "Oncorhynchus_kisutch")

NI.historic.plot <- ggplot(NI.sim.out.df.long) +
  facet_wrap(~initN) +
  geom_line(aes(x = years,
                y = mean,
                colour = factor(initN))) +
  geom_ribbon(aes(x = years,
                  ymin = lcl,
                  ymax = ucl,
                  fill = as.factor(initN)),
              alpha = 0.3) +
  geom_point(data = NIpops,
             aes(x = years,
                 y = popsize)) +
  
  scale_fill_manual(values = kispal) +
  scale_y_continuous(
    expand = expansion(mult = c(0)),
    limits = c(0, 10000),
    breaks = seq(0, 10000, 250)) +
  scale_x_continuous(limits = c(1987, 2017),
                     breaks = seq(1987, 2017, 1)) +
  
  theme_bw() +
  theme(legend.position = c(0.7, 0.15),
        legend.title = element_text(size = 8,
                                    colour = "black",
                                    face = "bold"),
        legend_text = element_text(size = 7),
        legend.direction = ("horizontal")) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(panel.grid = element_blank()) +
  
  theme(plot.title = element_blank()) +
  
  xlab("Years") +
  ylab("Population Size (breeding pairs") +
  labs(colour = expression("Pop Size 1987")) +
  theme(axis.title.x = element_blank()) + 
  theme(axis.ticks = element_line(size = 0.7)) +
  theme(axis.text.x = element_text(vjust = -0.5,
                                   size = 7,
                                   angle = 0)) +
  theme(axis.text.y = element_text(size = 7,
                                   margin = margin(r = 5))) +
  theme(axis.title.y = element_text(size = 10,
                                    face = "bold",
                                    margin = margin(r = 5))) +
  theme(axis.title.x = element_text(size = 10,
                                    face = "bold",
                                    margin = margin(t = 1,
                                                    r = 5,
                                                    b = 2),
                                    vjust = -1)) 
  
  




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



#----- Obtaining predictions that predict less than 750 curlew  -----
bettermatch <- which(NI.sim.out.3800[31 ,] < 750,
                     arr.ind = TRUE)
NI.sim.out.3800.match <- NI.sim.out.3800[, bettermatch]
NI.sim.out.3800.match.mean <- rowMeans(NI.sim.out.3800.match)
NI.sim.out.3800.match.lcl <- apply(NI.sim.out.3800.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.3800.match.ucl <- apply(NI.sim.out.3800.match, 1,
                                   quantile, probs = 0.9)

bettermatch <- which(NI.sim.out.4350[31 ,] < 750,
                     arr.ind = TRUE)
NI.sim.out.4350.match <- NI.sim.out.4350[, bettermatch]
NI.sim.out.4350.match.mean <- rowMeans(NI.sim.out.4350.match)
NI.sim.out.4350.match.lcl <- apply(NI.sim.out.4350.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.4350.match.ucl <- apply(NI.sim.out.4350.match, 1,
                                   quantile, probs = 0.9)


bettermatch <- which(NI.sim.out.4900[31 ,] < 750,
                     arr.ind = TRUE)
NI.sim.out.4900.match <- NI.sim.out.4900[, bettermatch]
NI.sim.out.4900.match.mean <- rowMeans(NI.sim.out.4900.match)
NI.sim.out.4900.match.lcl <- apply(NI.sim.out.4900.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.4900.match.ucl <- apply(NI.sim.out.4900.match, 1,
                                   quantile, probs = 0.9)

bettermatch <- which(NI.sim.out.5450[31 ,] < 750,
                     arr.ind = TRUE)
NI.sim.out.5450.match <- NI.sim.out.5450[, bettermatch]
NI.sim.out.5450.match.mean <- rowMeans(NI.sim.out.5450.match)
NI.sim.out.5450.match.lcl <- apply(NI.sim.out.5450.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.5450.match.ucl <- apply(NI.sim.out.5450.match, 1,
                                   quantile, probs = 0.9)



bettermatch <- which(NI.sim.out.6000[31 ,] < 750,
                     arr.ind = TRUE)
NI.sim.out.6000.match <- NI.sim.out.6000[, bettermatch]
NI.sim.out.6000.match.mean <- rowMeans(NI.sim.out.6000.match)
NI.sim.out.6000.match.lcl <- apply(NI.sim.out.6000.match, 1,
                                   quantile, probs = 0.1)
NI.sim.out.6000.match.ucl <- apply(NI.sim.out.6000.match, 1,
                                   quantile, probs = 0.9)


NI.sim.out.df.match <- data.frame(years = rep(years, times = 5),
                                  initN = rep(pop.1987, each = 31),
                                  mean = c(NI.sim.out.3800.match.mean,
                                           NI.sim.out.4350.match.mean,
                                           NI.sim.out.4900.match.mean,
                                           NI.sim.out.5450.match.mean,
                                           NI.sim.out.6000.match.mean),
                                  lcl = c(NI.sim.out.3800.match.lcl,
                                          NI.sim.out.4350.match.lcl,
                                          NI.sim.out.4900.match.lcl,
                                          NI.sim.out.5450.match.lcl,
                                          NI.sim.out.6000.match.lcl),
                                  ucl = c(NI.sim.out.3800.match.ucl,
                                          NI.sim.out.4350.match.ucl,
                                          NI.sim.out.4900.match.ucl,
                                          NI.sim.out.5450.match.ucl,
                                          NI.sim.out.6000.match.ucl)
)

write.table(NI.sim.out.df.match,
            "data/NI_historic_1987_match750_1.csv",
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
        "data/ROI_historic1985_sims.RDS")
# version 2 of ROI 1985 base saved with breed.fail 0.5, ROI.K = 25000



#-----  ROI Plotting  -----
ROI.sim.out <- readRDS("data/ROI_historic1985_sims.RDS")

roi.names <- rep(as.character(ROI.pop.1985), (4725/5))
names(ROI.sim.out.hist) <- roi.names

years.roi.hist <- seq(1985, 2015, 1)

ROI.sim.out.6000 <- ROI.sim.out.hist[names(ROI.sim.out.hist) == "6000"]
ROI.sim.out.6000 <- do.call(cbind, ROI.sim.out.6000)
ROI.sim.out.6000.mean <- rowMeans(ROI.sim.out.6000)
ROI.sim.out.6000.lcl  <- apply(ROI.sim.out.6000, 1, quantile, probs = 0.1)
ROI.sim.out.6000.ucl  <- apply(ROI.sim.out.6000, 1, quantile, probs = 0.9)

ROI.sim.out.6500 <- ROI.sim.out.hist[names(ROI.sim.out.hist) == "6500"]
ROI.sim.out.6500 <- do.call(cbind, ROI.sim.out.6500)
ROI.sim.out.6500.mean <- rowMeans(ROI.sim.out.6500)
ROI.sim.out.6500.lcl  <- apply(ROI.sim.out.6500, 1, quantile, probs = 0.1)
ROI.sim.out.6500.ucl  <- apply(ROI.sim.out.6500, 1, quantile, probs = 0.9)

ROI.sim.out.7000 <- ROI.sim.out.hist[names(ROI.sim.out.hist) == "7000"]
ROI.sim.out.7000 <- do.call(cbind, ROI.sim.out.7000)
ROI.sim.out.7000.mean <- rowMeans(ROI.sim.out.7000)
ROI.sim.out.7000.lcl  <- apply(ROI.sim.out.7000, 1, quantile, probs = 0.1)
ROI.sim.out.7000.ucl  <- apply(ROI.sim.out.7000, 1, quantile, probs = 0.9)

ROI.sim.out.7500 <- ROI.sim.out.hist[names(ROI.sim.out.hist) == "7500"]
ROI.sim.out.7500 <- do.call(cbind, ROI.sim.out.7500)
ROI.sim.out.7500.mean <- rowMeans(ROI.sim.out.7500)
ROI.sim.out.7500.lcl  <- apply(ROI.sim.out.7500, 1, quantile, probs = 0.1)
ROI.sim.out.7500.ucl  <- apply(ROI.sim.out.7500, 1, quantile, probs = 0.9)

ROI.sim.out.8000 <- ROI.sim.out.hist[names(ROI.sim.out.hist) == "8000"]
ROI.sim.out.8000 <- do.call(cbind, ROI.sim.out.8000)
ROI.sim.out.8000.mean <- rowMeans(ROI.sim.out.8000)
ROI.sim.out.8000.lcl  <- apply(ROI.sim.out.8000, 1, quantile, probs = 0.1)
ROI.sim.out.8000.ucl  <- apply(ROI.sim.out.8000, 1, quantile, probs = 0.9)


# Combine into tibble


ROI.sim.out.df.long <- data.frame(years = rep(years, times = 5),
                             initN = rep(ROI.pop.1985, each = 31),
                             mean = c(ROI.sim.out.6000.mean,
                                      ROI.sim.out.6500.mean,
                                      ROI.sim.out.7000.mean,
                                      ROI.sim.out.7500.mean,
                                      ROI.sim.out.8000.mean),
                             lcl = c(ROI.sim.out.6000.lcl,
                                     ROI.sim.out.6500.lcl,
                                     ROI.sim.out.7000.lcl,
                                     ROI.sim.out.7500.lcl,
                                     ROI.sim.out.8000.lcl),
                             ucl = c(ROI.sim.out.6000.ucl,
                                     ROI.sim.out.6500.ucl,
                                     ROI.sim.out.7000.ucl,
                                     ROI.sim.out.7500.ucl,
                                     ROI.sim.out.8000.ucl)
)


write.table(ROI.sim.out.df.long,
            "data/ROI_sim_out_historic1985_means_1.csv",
            sep = ",",
            row.names = FALSE)

# divide by 2 to get breeding pairs



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




#----- Obtaining predictions that predict less than 400 curlew  -----
bettermatch.6000 <- which(ROI.sim.out.6000[31 ,] < 400,
                          arr.ind = TRUE)
ROI.sim.out.6000.match <- ROI.sim.out.6000[, bettermatch.6000]
ROI.sim.out.6000.match.mean <- rowMeans(ROI.sim.out.6000.match)
ROI.sim.out.6000.match.lcl <- apply(ROI.sim.out.6000.match, 1,
                                    quantile, probs = 0.1)
ROI.sim.out.6000.match.ucl <- apply(ROI.sim.out.6000.match, 1,
                                    quantile, probs = 0.9)

bettermatch.6500 <- which(ROI.sim.out.6500[31 ,] < 400,
                          arr.ind = TRUE)
ROI.sim.out.6500.match <- ROI.sim.out.6500[, bettermatch.6500]
ROI.sim.out.6500.match.mean <- rowMeans(ROI.sim.out.6500.match)
ROI.sim.out.6500.match.lcl <- apply(ROI.sim.out.6500.match, 1,
                                    quantile, probs = 0.1)
ROI.sim.out.6500.match.ucl <- apply(ROI.sim.out.6500.match, 1,
                                    quantile, probs = 0.9)

bettermatch.7000 <- which(ROI.sim.out.7000[31 ,] < 400,
                          arr.ind = TRUE)
ROI.sim.out.7000.match <- ROI.sim.out.7000[, bettermatch.7000]
ROI.sim.out.7000.match.mean <- rowMeans(ROI.sim.out.7000.match)
ROI.sim.out.7000.match.lcl <- apply(ROI.sim.out.7000.match, 1,
                                    quantile, probs = 0.1)
ROI.sim.out.7000.match.ucl <- apply(ROI.sim.out.7000.match, 1,
                                    quantile, probs = 0.9)

bettermatch.7500 <- which(ROI.sim.out.7500[31 ,] < 400,
                          arr.ind = TRUE)
ROI.sim.out.7500.match <- ROI.sim.out.7500[, bettermatch.7500]
ROI.sim.out.7500.match.mean <- rowMeans(ROI.sim.out.7500.match)
ROI.sim.out.7500.match.lcl <- apply(ROI.sim.out.7500.match, 1,
                                    quantile, probs = 0.1)
ROI.sim.out.7500.match.ucl <- apply(ROI.sim.out.7500.match, 1,
                                    quantile, probs = 0.9)

bettermatch.8000 <- which(ROI.sim.out.8000[31 ,] < 400,
                          arr.ind = TRUE)
ROI.sim.out.8000.match <- ROI.sim.out.8000[, bettermatch.8000]
ROI.sim.out.8000.match.mean <- rowMeans(ROI.sim.out.8000.match)
ROI.sim.out.8000.match.lcl <- apply(ROI.sim.out.8000.match, 1,
                                    quantile, probs = 0.1)
ROI.sim.out.8000.match.ucl <- apply(ROI.sim.out.8000.match, 1,
                                    quantile, probs = 0.9)




ROI.sim.out.df.match <- data.frame(years = rep(years, times = 5),
                                  initN = rep(ROI.pop.1985, each = 31),
                                  mean = c(ROI.sim.out.6000.match.mean,
                                           ROI.sim.out.6500.match.mean,
                                           ROI.sim.out.7000.match.mean,
                                           ROI.sim.out.7500.match.mean,
                                           ROI.sim.out.8000.match.mean),
                                  lcl = c(ROI.sim.out.6000.match.lcl,
                                          ROI.sim.out.6500.match.lcl,
                                          ROI.sim.out.7000.match.lcl,
                                          ROI.sim.out.7500.match.lcl,
                                          ROI.sim.out.8000.match.lcl),
                                  ucl = c(ROI.sim.out.6000.match.ucl,
                                          ROI.sim.out.6500.match.ucl,
                                          ROI.sim.out.7000.match.ucl,
                                          ROI.sim.out.7500.match.ucl,
                                          ROI.sim.out.8000.match.ucl)
)

write.table(ROI.sim.out.df.match,
            "data/ROI_sims_1985_match400.csv",
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



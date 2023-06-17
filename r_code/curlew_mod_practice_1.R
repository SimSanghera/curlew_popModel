######  Simple Curlew Population Viability Analysis #####

# Model adapted from Shoemaker and Oppel


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
NI.pop.sizes    <- seq(250, 560, 10)      # NI pop range from 250 - 500 
                                          # with 10 pair intervals

ROI.pop.sizes   <- seq(100, 200, 10)      # ROI pop range

I.popsizes      <- seq(100, 560, 20)      # all ireland pop size

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

intro.years     <- seq(3, 10, 1)            # the number of years in which
                                          # birds will be introduced


# Specify parameters for Pop Viability Analysis
n.years   <- 30                           # number of years to simulate
n.reps    <- 1                          # number of simulations
K         <- 5000                         # carrying capacity of site
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

extinct.cutoff  <- c(10, 20, 50, 100)     # below this number, pop is extinct
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
  S.jw, S.ad
)

bird.vr <- list(prod.rate = mean(prod.rate),
                S.ad  = mean(S.ad),
                S.jw  = mean(S.jw),
                sx.r = sex.ratio)

##########################################################################
bird.vr.range <- list(prod.rate,
                      S.ad,
                      S.jw,
                      sx.r = sex.ratio)

les.mat <- matrix(sapply(bird.matrix,
                         eval,
                         bird.vr,
                         NULL),
                  nrow = sqrt(length(bird.matrix)),
                  byrow = TRUE)

names(bird.vr.range) <- c("prod.rate",
                          "S.ad",
                          "S.jw",
                          "sx.r")

# won't work as matrix can't hold multiple values in each cell, so need
# to repeat simulations and run through full set of combinations
# individually
les.mat.range <- matrix(sapply(bird.matrix,
                               eval,
                               bird.vr.range,
                               NULL),
                        nrow = sqrt(length(bird.matrix)),
                        byrow = TRUE)

# projections have 50 wild juvs, 10 released juvs, 100 adults
# so we will again need the number of released juveniles to be a range from
# 0 - 50?
projections <- pop.projection(les.mat,
                              n = c(50, 10, 100),     
                              iterations = 50)
projections$lambda

projections.range <- pop.projection(les.mat.range,
                                    n = c(50, 0, 250),
                                    iterations = 50)
############################################################################


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


# Stochastic PVA

R.max <- 1                        # maximum pop growth rate

PVA <- function(
  n.reps,
  n.years,
  init.N,
  K,
  catastrophe.p,
  catastrophe.i,
  F.loop, # obtained from simul.in demographics of fecundity
  S.jw.loop,
  S.ad.loop,
  breed.fail
#  sex.ratio  - remember to add , back in after breed.fail
) {
  
  # Create array to store population values
  # rows = years, cols = reps
  Pop.Array <- array(0,
                     dim = c((n.years + 1),
                             n.reps))
  
  # Start loop of reps
  for (rep in 1:n.reps) {
    
    # initial abundance minus number of head-start birds
    # but we will have a seq of head-starting scenarios, so won't need
    # to select from a poisson distribution
    Pop.Array[1, rep]   <- init.N*2 #- rpois(1, capt.fail * init.N)
    
    # Loop through the years
    for (y in 2:(n.years + 1)) {
      
      # Create matrix with vital rates
      # Is fecundity population size dependent? 
      # F.Year <- F.loop * ifelse(Pop.Array[y - 1, rep] < 20, 0.25, 1)
      
      # Ruined breeding season due to weather/disturbance/etc...
      # Possibly count number of poor breeding years in survey data and
      #   obtain a more accurate probability fo this occurring
      # f.loop comes from simul.in table
      
      F.Year <- ifelse(rbinom(1, 1, breed.fail) == 1,
                       F.loop * 0.5, 
                       F.loop)
      
      bird.vr <- list(prod.rates = F.Year,
                      S.ad = S.ad.loop,
                      S.jw = S.jw.loop #,
                      # sx.r = sex.ratio
                      )
      
      # les.mat <- matrix(sapply(bird.matrix,
      #                          eval,
      #                          bird.vr,
      #                          NULL),
      #                   nrow = sqrt(length(bird.matrix)),
      #                   byrow = TRUE)
      a <- matrix(sapply(bird.matrix,
                         eval,
                         bird.vr,
                         NULL),
                  nrow = sqrt(length(bird.matrix)),
                  byrow = TRUE)
      
      pop.size  <- c(((init.N*2)/3), ((init.N*2)/3)*2) # why divide by 3?
                        # had to multiply init.N by 2 to match above code
      
      projections <- pop.projection(a,
                                    #les.mat
                                    n = pop.size,
                                    iterations = 30)
      
      # Stochastic Ricker Model
      R.max <- projections$lambda      # max growth rate (max lambda)
      
      ## stochasticity & density dependence
      # calculate abundance based on Ricker model
      # - rounded to integer, set to min of 0
      next.year <- max(0, 
                       trunc(Ricker(Pop.Array[y-1, rep])))
      
      # Catastrophe - check prob and then calculate impact on pop
      if(runif(1) < catastrophe.p) next.year <- next.year * catastrophe.i
      
      # catastrophe when in release program
      # if(y < intro.years) next.year <- next.year + init.N
      
      # Obtain year abundance
      Pop.Array[y, rep] <- next.year
      
    }
    
  }
  
  # Return
  return(Pop.Array)
  
}


#-----  Calculate Proportion of Sims where Species Goe Extinct  -----
extinction.bysim  <- function(simdata,
                              threshold) {
  sum(apply(default,
            2,
            function(t) min(t) < threshold)) / ncol(simdata)
  ## extinction is defined as < threshold number of birds
  
}


#-----  Loop over Combinations of Demography  -----

# Table of all combinations of demographic parameters


simul.in <- expand.grid(I.popsizes,
                        S.ad,
                        S.jw,
                        prod.rates)  # for release, add intro years

dim(simul.in)
names(simul.in) <- c("pop.size",
                     "S.ad",
                     "S.jw",
                     "Productivity Rate (F)")

sim.out <- data.frame()
demographics <- list()    # to store annual demographics
                          # total pop
                          # age class numbers, prod rate, survival, 
                          # number of 1st year not returning
                          # number of released birds

# Boost speed by setting up parallel processors
cl <- makeCluster(4)
registerDoParallel(cl, cores = 4)


#-----  Start parallel loop -----

sim.out <- foreach(s = c(1:dim(simul.in)[1]),
                   .packages = "popbio",
                   .combine = rbind) %dopar% {
                     
                     init.N <- simul.in[s, 1]
                     default <- PVA(n.reps,
                                    n.years,
                                    init.N,
                                    K,
                                    catastrophe.p,
                                    catastrophe.i,
                                    F.loop = simul.in[s, 4],
                                    S.jw.loop = simul.in[s, 3],
                                    S.ad.loop = simul.in[s, 2],
                                    breed.fail)
                     
                     # Calculating mean population growth rate
                     out <- simul.in[s, ]
                     bird.vr <- list(prod.rates = simul.in[s, 4],
                                     S.ad = simul.in[s, 2],
                                     S.jw = simul.in[s, 3])
                     
                     a <- matrix(sapply(bird.matrix,
                                        eval,
                                        bird.vr,
                                        NULL),
                                 nrow = sqrt(length(bird.matrix)),
                                 byrow = TRUE)
                     pop.size <- c(0, (init.N/2))
                     projections <- pop.projection(a, 
                                                   n = pop.size,
                                                   iterations = 30)
                     out$lambda <- projections$lambda
                     
                     
                     # Calculate Extinction Probability
                     final.out <- data.frame()
                     for (t in extinct.cutoff) {
                       out$ExtThresh <- t
                       out$outcome <- extinction.bysim(default,
                                                           t)
                       final.out <- rbind(final.out, out)
                       
                     }
                     
                     return(final.out)
                     
                   }

# stop parallel processing
stopCluster(cl)

write.table(sim.out,
            "/model_output/curlew_base_model_trial_1.csv",
            sep = ",",
            row.names = FALSE)



#-----  Summarise Output  -----

heed(sim.out)
hist(sim.out$outcome)


# Plot simulation output
sim.out.2 <- sim.out %>%
  group_by(pop.size,
           ExtThresh) %>%
  summarise(prop.ext = mean(outcome),
            lcl = quantile(outcome, 0.1),
            ucl = quantile(outcome, 0.9)) %>%
  mutate(extinct = ifelse(ExtThresh == 1,
                          "< 1 bird",
                          paste("<", ExtThresh, "birds")))

ext.risk.plot <- ggplot(sim.out.2) + 
  geom_line(aes(x = pop.size,
                y = prop.ext,
                color = extinct)) +
  geom_ribbon(aes(x = pop.size,
                  ymin = lcl,
                  ymax = ucl,
                  fill = extinct),
              alpha = 0.3) +
  xlab("Number of birds introduced per year") +
  ylab("Probability of curlew extinct within x years") +
  labs(color = expression(paste("Population n extinct when")),
       fill = expression(paste("Population n extinct when"))) +
  theme(panel.background=element_rect(fill="white", colour="black"),
        legend.background = element_rect(),
        legend.title = element_text(size=16),
        legend.text=element_text(size=14),
        legend.position=c(0.8,0.82),
        axis.text=element_text(size=16, color="black"), 
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=16, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())
ext.risk.plot

ggsave("figures/curlew_PVA_trial_extprob.png", width = 9, height = 6)


  
  




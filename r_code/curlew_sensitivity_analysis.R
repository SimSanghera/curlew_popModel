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
F.Year           = simul.in.simple[s, 7]
S.ad.loop        = simul.in.simple[s, 2]
S.jw.loop        = simul.in.simple[s, 3]
S.jhs.loop       = simul.in.simple[s, 4]
S.imm.loop       = simul.in.simple[s, 5]
S.sub.loop       = simul.in.simple[s, 6]
philopatry       = simul.in.simple[s, 8]
sitefidelity     = simul.in.simple[s, 9]
hs.numbers       = simul.in.simple[s, 10]
hs.years         = simul.in.simple[s, 11]
pop.size <- c(0, 0, 0, 0, (init.N*2))



# Changes to philopatry
#  3 values - 0.3, 0.5, 0.7

# Changes to survival, as noted above

# Changes to return, as noted above

# F.Year just keep as prod rates, change as above

# Constants: initN, hs.numbers, hs.years



# Process: get down to projections part of code, then run sens and elasticities
#  Extract data and put into table
# Then re-run with new param values

#----- Philopatry -----
# 0.3
w <- eigen(a)$vectors
v <- Conj(solve(w))
senmat <- Re(v[1, ] %*% t(w[, 1]))
emat <- (1/(Re(eigen(a)$values[1]))) * senmat * a

philopatry = 0.3

phil_low_sens <- sensitivity(a, zero = FALSE)
phil_low_elas <- elasticity(a)

phil.p.low.projections <- pop.projection(a,
                                          n = pop.size,
                                          iterations  = 30)

phil.low.eig <- eigen.analysis(a)

# 0.5
philopatry = 0.5
phil_mid_sens <- sensitivity(a, zero = FALSE)
phil_mid_elas <- elasticity(a)

phil.p.med.projections <- pop.projection(a,
                                          n = pop.size,
                                          iterations  = 30)

phil.med.eig <- eigen.analysis(a)


# 0.7
philopatry = 0.7
phil_high_sens <- sensitivity(a, zero = FALSE)
phil_high_elas <- elasticity(a)

phil.p.high.projections <- pop.projection(a,
                                          n = pop.size,
                                          iterations  = 30)

phil.high.eig <- eigen.analysis(a)
image(phil.high.eig$elasticities)


philopatry = 1
phil.massive.sens <- sensitivity(a, zero = FALSE)
phil.massive.elas <- elasticity(a)

phil.p.massive.projections <- pop.projection(a,
                                          n = pop.size,
                                          iterations  = 30)

phil.massive.eig <- eigen.analysis(a)
image(phil.high.eig$elasticities)



#-----  Adult Survival  -----
# keep phil.p at 0.5
# 0.78 
S.ad.loop <- 0.78
s.ad.low.sens <- sensitivity(a)
s.ad.low.elas <- elasticity(a)

s.ad.low.proj <- pop.projection(a,
                                n = pop.size,
                                iterations = 30)

s.ad.low.eig <- eigen.analysis(a)

# 0.84
S.ad.loop <- 0.84
s.ad.med.sens <- sensitivity(a)
s.ad.mid.elas <- elasticity(a)

s.ad.med.proj <- pop.projection(a,
                                n = pop.size,
                                iterations = 30)

s.ad.med.eig <- eigen.analysis(a)


# 0.86
s.ad.hm.elas <- elasticity(a)


# 0.90
S.ad.loop <- 0.90
s.ad.high.sens <- sensitivity(a)
s.ad.high.elas <- elasticity(a)

s.ad.high.proj <- pop.projection(a,
                                n = pop.size,
                                iterations = 30)

s.ad.high.eig <- eigen.analysis(a)


#----- JW Fledge Survival  -----
# keep phil.p at 0.5, S.ad = 0.78
# 0.2
S.jw.loop <- 0.2
s.jw.low.sens <- sensitivity(a)
s.jw.low.elas <- elasticity(a)

s.jw.low.proj <- pop.projection(a,
                                n = pop.size,
                                iterations = 30)

s.jw.low.eig <- eigen.analysis(a)

# 0.32
S.jw.loop <- 0.32
s.jw.med.sens <- sensitivity(a)
s.jw.mid.elas <- elasticity(a)

s.jw.med.proj <- pop.projection(a,
                                n = pop.size,
                                iterations = 30)

s.jw.med.eig <- eigen.analysis(a)


# 0.37
S.jw.loop <- 0.37
s.jw.high.sens <- sensitivity(a)
s.jw.high.elas <- elasticity(a)

s.jw.high.proj <- pop.projection(a,
                                 n = pop.size,
                                 iterations = 30)

s.jw.high.eig <- eigen.analysis(a)

S.jw.loop <- 0.87
s.jw.massive.sens <- sensitivity(a)
s.jw.massive.elas <- elasticity(a)

s.jw.massive.proj <- pop.projection(a,
                                 n = pop.size,
                                 iterations = 30)

s.jw.massive.eig <- eigen.analysis(a)




#----- Imm Survival  -----
# keep phil.p at 0.5, S.ad = 0.78, s.jw = 0.2
# 0.4
S.imm.loop <- 0.4
s.imm.low.sens <- sensitivity(a)
s.imm.low.elas <- elasticity(a)

s.imm.low.proj <- pop.projection(a,
                                n = pop.size,
                                iterations = 30)

s.imm.low.eig <- eigen.analysis(a)

# 0.5
S.imm.loop <- 0.5
s.imm.med.sens <- sensitivity(a)
s.imm.mid.elas <- elasticity(a)

s.imm.med.proj <- pop.projection(a,
                                n = pop.size,
                                iterations = 30)

s.imm.med.eig <- eigen.analysis(a)


# 0.7
S.imm.loop <- 0.7
s.imm.high.sens <- sensitivity(a)
s.imm.high.elas <- elasticity(a)

s.imm.high.proj <- pop.projection(a,
                                 n = pop.size,
                                 iterations = 30)

s.imm.high.eig <- eigen.analysis(a)



#----- Sub Survival  -----
# keep phil.p at 0.5, S.ad = 0.78
# 0.78
S.sub.loop <- 0.78
s.sub.low.sens <- sensitivity(a)
s.sub.low.elas <- elasticity(a)

s.sub.low.proj <- pop.projection(a,
                                n = pop.size,
                                iterations = 30)

s.sub.low.eig <- eigen.analysis(a)

# 0.84
S.sub.loop <- 0.84
s.sub.med.sens <- sensitivity(a)
s.sub.mid.elas <- elasticity(a)

s.sub.med.proj <- pop.projection(a,
                                n = pop.size,
                                iterations = 30)

s.sub.med.eig <- eigen.analysis(a)


# 0.9
S.sub.loop <- 0.90
s.sub.high.sens <- sensitivity(a)
s.sub.high.elas <- elasticity(a)

s.sub.high.proj <- pop.projection(a,
                                 n = pop.size,
                                 iterations = 30)

s.sub.high.eig <- eigen.analysis(a)



#----- Prod Rates  -----
# keep phil.p at 0.5, S.ad = 0.78
# 0.78
F.Year <- 0.14
prod.low.sens <- sensitivity(a)
prod.low.elas <- elasticity(a)

prod.low.proj <- pop.projection(a,
                                 n = pop.size,
                                 iterations = 30)

prod.low.eig <- eigen.analysis(a)

# 0.84
F.Year<- 0.32
prod.med.sens <- sensitivity(a)
prod.mid.elas <- elasticity(a)

prod.med.proj <- pop.projection(a,
                                 n = pop.size,
                                 iterations = 30)

prod.med.eig <- eigen.analysis(a)


# 0.9
F.Year <- 0.5
prod.high.sens <- sensitivity(a)
prod.high.elas <- elasticity(a)

prod.high.proj <- pop.projection(a,
                                 n = pop.size,
                                 iterations = 30)

prod.high.eig <- eigen.analysis(a)



#----- initN  -----
# keep phil.p at 0.5, S.ad = 0.78
# 0.78
init.N <- 100
N.low.sens <- sensitivity(a)
N.low.elas <- elasticity(a)

N.low.proj <- pop.projection(a,
                                n = pop.size,
                                iterations = 30)

N.low.eig <- eigen.analysis(a)

# 150
init.N <- 150
N.med.sens <- sensitivity(a)
N.mid.elas <- elasticity(a)

N.med.proj <- pop.projection(a,
                                n = pop.size,
                                iterations = 30)

N.med.eig <- eigen.analysis(a)


# 200
init.N <- 200
N.high.sens <- sensitivity(a)
N.high.elas <- elasticity(a)

N.high.proj <- pop.projection(a,
                                 n = pop.size,
                                 iterations = 30)

N.high.eig <- eigen.analysis(a)



#----- Return Rate Adults -----
# 0.75

sitefidelity = 0.75

return.low.sens <- sensitivity(a, zero = FALSE)
return.low.elas <- elasticity(a)

return.low.projections <- pop.projection(a,
                                         n = pop.size,
                                         iterations  = 30)

return.low.eig <- eigen.analysis(a)

# 0.85
sitefidelity = 0.85
return.med.sens <- sensitivity(a, zero = FALSE)
return.mid.elas <- elasticity(a)

return.med.projections <- pop.projection(a,
                                         n = pop.size,
                                         iterations  = 30)

return.med.eig <- eigen.analysis(a)


# 0.95
sitefidelity = 0.95
return.high.sens <- sensitivity(a, zero = FALSE)
return.high.elas <- elasticity(a)

return.high.projections <- pop.projection(a,
                                          n = pop.size,
                                          iterations  = 30)

return.high.eig <- eigen.analysis(a)



#----- Combinations -----
# Combine imm and sub survival & phil
S.imm.loop <- 0.9
S.sub.loop <- 0.9
S.jw.loop <- 0.2
phil.p = 1

combo.imm.sub.phil.elas <- elasticity(a)
combo.imm.sub.phil.proj <- pop.projection(a,
                                           n = pop.size,
                                           iterations = 30)

combo.imm.sub.phil.eig <- eigen.analysis(a)




#----- Plot -----
# Extract Lambda values from each and plot against low med high for each
#   param


# Params first
Philopatry        <- c(0.3, 0.5, 0.7)
AdultSurvival     <- c(0.7, 0.84, 0.90)
FledgeSurvival    <- c(0.2, 0.32, 0.37)
ImmSurvival       <- c(0.4, 0.5, 0.7)
SubSurvival       <- c(0.78, 0.84, 0.90)
ProdRates         <- c(0.14, 0.32, 0.5)
initN             <- c(100, 150, 200)
AdultReturn       <- c(0.75, 0.85, 0.95)


# Lambdas
PhilLambda         = c(0.7438743, 0.7457557, 0.7476104)
S.adLambda         = c(0.7457557, 0.808257, 0.8581207)
S.jwLambda         = c(0.7457557, 0.748528, 0.7496661)
S.immLambda        = c(0.7457557, 0.7469179, 0.7492117)
S.subLambda        = c(0.7457557, 0.7461142, 0.7464716)
ProdRatesLambda    = c(0.7457557, 0.7516267, 0.7572506)
initNLambda        = c(0.7457557, 0.7457557, 0.7457557)
ReturnLambda       = c(0.5909592, 0.6682914, 0.7457557)

Lambda <- c(PhilLambda, S.adLambda,
            S.jwLambda, S.immLambda,
            S.subLambda, ProdRatesLambda,
            initNLambda, ReturnLambda)

Perc_Diff <- ((Lambda - 0.7457557)/0.7457557) * 100


values <- c("low", "med", "high")
parameters <- c("Philopatry",
                "Adult Survival",
                "Fledgling Survival",
                "Immature Survival",
                "Subadult Survival",
                "Prod Rates",
                "Initial Pop Size",
                "Adult Return Rate")

elas.df <- tibble(
  Parameters = rep(parameters, each = 3),
  Values     = rep(values, times = 8),
  Lambda     = Lambda,
  PercDiff  = Perc_Diff
)

write.table(elas.df, 
            "data/elasticity_analysis.csv",
            sep = ",",
            row.names = FALSE)


# Plot Lambdas
library(wesanderson)
foxPal <- wes_palette("FantasticFox1")



elas.plot <- ggplot(data = elas.df,
                    aes(x = Parameters,
                        y = Lambda,
                        group = Values,
                        fill = Values)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  geom_text(aes(group = Values,
                label = round(PercDiff, digits = 2)),
            position = position_dodge(.9),
            vjust = -1) +
  scale_fill_manual(values = foxPal) +
  geom_hline(yintercept = 0.7457557) +
  xlab("Parameters") +
  ylab("Population Growth Rate\n") +
  scale_y_continuous(expand = expansion(mult = c(0)),
                     limits = c(0.0, 1.0),
                     breaks = seq(0.0, 1.0, 0.1)) +
  coord_cartesian(ylim = c(0.0, 1.0)) + 
  labs(colour = "Parameters") +
  theme_bw() +
  theme(axis.title = element_text(size = 32),
        axis.title.x = element_text(vjust = -0.5),
        axis.text = element_text(size = 16, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        legend.title = element_text(size = 20),
        legend.position = c(0.95, 0.91),
        legend.justification = "right",
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size = 20))



###################
stage.vector.plot(projA$stage.vectors)


plot(projections$pop.sizes,
     pch = 16,
     type = "b",
     xlab = "Year",
     ylab = "Population Size")


# | One of the most powerful aspects of using population matrix models is that we can use the mathematical techniques of linear
# | algebra to learn about the dynamics of the population we are modeling and estimate some informative demographic parameters. We
# | can get an estimate of lambda that is independent of the length of time that we model population growth (since lambda should
#                                                                                                            | stabilize as long as all our age-specific survival and fertility probabilities are constant). We can determine the stable
# | age/stage distribution. And we can see how different elements of our matrix (the vital rate) affect our growth rate lambda.
# 
# | The function eigen.analysis() can accomplish all these goals. It takes only one argument - our matrix. If we called our matrix
# | A, we would simply use the command eigen.analysis(A). We will give the results of this analysis a name - eigA. Use this code to
# | perform the eigen.analysis - eigA=eigen.analysis(A)

eigA = eigen.analysis(a)
eigA



# | The output of eigen.analysis has a lot of useful information. First we see by $lambda1 the value of lambda the population will
# | achieve when it is stabilized. $stable.stage shows the stable stage/age distribution. Then we see two matrices of values - one
# | for sensitivities and one for elasticities. These both indicate the importance of each element of the matrix (the fertility
# | rates, survival probabilities, etc.) in determining the growth rate of the population, lambda. Technically, sensitivities
# | indicate how a small absolute change in a particular age-specific fertility value (Fi) or survival probability (Pi) will affect
# | lambda when all other matrix elements are held constant. Elasticities are an even more useful measure - they indicate how a
# | proportional change in an Fi or Pi value will affect lambda when all other matrix elements are held constant. Because
# | elasticities sum to 1, they are easier to compare and interpret than sensitivities.Finally, the $repro.value (reproductive
# | value) indicates the expected contribution of an individual of each age or stage class to both current and future reproduction.
# | These reproductive values are scaled relative to the values of the first age class, which is set to 1.


image2(eigA$elasticities)




####################
giraffe <- matrix(c(0, 0, 0.24, 0.57, 0, 0, 0, 0.79, 0.84), nrow=3, byrow=TRUE,
                  dimnames=list(c("calf","subadult","adult"),
                                c("calf","subadult","adult")))

mat <- giraffe

mat

sensitivity(mat, zero = FALSE)

elasticity(mat)


# Create values for transition rates
vals <- c(1.01, 1.05, 1.1)

# Create list of results for each transition probability
results <- list()
for (x in 1:length(vals)) {
  testlam <- matrix(0, nrow=nrow(mat), ncol=nrow(mat))
  
  for (i in c(1:nrow(mat))) {
    for (j in c(1:nrow(mat))) {
      if (mat[i,j] == 0) {testlam[i,j] <- 0} else {
        tempmat <- mat
        tempmat[i,j] <- mat[i,j]*vals[x]
        testlam[i,j] <- Re(eigen(tempmat)$values[1])
      }
    }
  }
  results[[x]] <- testlam
}


output <- matrix(NA, nrow=length(vals), ncol=length(which(as.vector(t(results[[1]])) >0)))
for (w in c(1: length(vals))) {
  output[w,] <- as.vector(t(results[[w]]))[which(as.vector(t(results[[w]])) >0)]
}


lower <- 0.6

par(mar=c(5, 5, 2, 2))
barplot(output-lower, beside=TRUE, ylim=c(lower, 1.12), ylab="Population growth rate",
        col=c("#007765","#6AA341","gray90"), offset= lower)
abline(h=lower, lwd=2)
abline(h=1.0, lwd=2, col="blue")
abline(h= (Re(eigen(a)$values[1])), lwd=2, col="orange")
legend("topleft", inset=c(0.01,0.01), c("+1%","+5%","+10%"),
       fill=c("#007765","#6AA341","gray90"), cex=0.8)
text(x=c(3.7,7.7,11.65,15.8, 19.6), y=0.78,
     labels=c("AdultProd","FledgeSur","ImmSurv", "SubSurv", "AdultSurv"),
     pos=2, xpd=TRUE, cex=0.9)








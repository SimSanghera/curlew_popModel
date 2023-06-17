#-----  Notes -----
# phil.p = philopatry probability. It turns out curlew may not be as philopatric
#     as everyone assumed and there is no data on this. This would multiply
#     the survival of the S.jw and S.jhs birds. You can see in the stage
#     matrix below.
# return.p = return to site probability. This is typically high in adult 
#     curlew, once they have bred. So i have set that between 0.75-0.95.
#     I could reduce this to just 0.7, 0.8, 0.9. I don't think the extra .05
#     will make too much difference.

# I think the big ones are productivity rate (prod.rates) and numbers of 
#   introduced HS birds (hs.pop). Productivity is very low in curlew, so
#   jumping too much of the low end might over-estimate their survival success.
#   

# Steffen suggested subsetting the data - so maybe run scenarios with 
#   high levels of return and low levels?
# Same for philopatry


# Also split into several base models - with low levels of philopatry to high
#   4 models - 0.3 - 0.5, 0.7 - 0.9
# Low and high levels of return - 3 models
# 


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


# Population Size
NI.pop.size   <- seq(250, 500, 20)      # breeding pairs. Closer to 250.

ROI.pop.size  <- seq(100, 200, 20)


# Original Survival Rates & Curlew Related Parameters
S.ad          <- seq(0.68, 0.92, 0.02)  # Cook et al., 2021
                                        # Roodbergen et al., 2021

S.jw          <- seq(0.20, 0.40, 0.02)  # Wild bird survival
                                        # Cook et al., 2021

S.jhs         <- seq(0.20, 0.40, 0.02)  # HS bird survival
                                        # Assume same as wild for now
                                        # Geoff Hilton (pers comms)

S.sub         <- seq(0.68, 0.92, 0.04)  # subadult survival

S.imm         <- seq(0.20, 0.50, 0.05)  # survival of immature (age fledge - 1)

prod.rates    <- seq(0.14, 0.52, 0.02)  # Grant et al., 1999
                                        # Harris et al., 2019
                                        # Zielonka et al., 2019
                                        # Irish Breeding Bird Report

hs.pop        <- seq(6, 60, 3)          # this number represents the number of
                                        # individual birds added, but split
                                        # by clutch
                                        # assume 3 eggs per clutch
                                        # i.e. 2 clutches to 20

capt.fail     <- 0.2                    # The number of HS birds that
                                        # fail during rearing & can't be
                                        # released

intro.years   <- seq(1, 5, 1)           # years in which hs birds will be
                                        # introduced

phil.p        <- seq(0.30, 0.90, 0.1)   # probability of first years return
                                        # to release site

return.p      <- seq(0.75, 0.95, 0.05)  # probability of breeding birds 
                                        # returning to release site




#-----  Reduced Values -----  
NI.pop.size   <- seq(200, 300, 25)

ROI.pop.size  <- seq(100, 150, 25)


# Reducing range of Survival Rates & Curlew Related Parameters
S.ad          <- seq(0.78, 0.92, 0.04)  # Cook et al., 2021
                                        # Roodbergen et al., 2021

S.jw          <- seq(0.25, 0.40, 0.05)  # Wild bird survival
                                        # Cook et al., 2021

S.jhs         <- seq(0.25, 0.40, 0.05)  # HS bird survival
                                        # Assume same as wild for now
                                        # Geoff Hilton (pers comms)

S.sub         <- seq(0.78, 0.92, 0.04)  # subadult survival

S.imm         <- seq(0.40, 0.70, 0.1)   # survival of immature (age 1 - 2)

prod.rates    <- seq(0.14, 0.52, 0.06)  # Grant et al., 1999
                                        # Harris et al., 2019
                                        # Zielonka et al., 2019
                                        # Irish Breeding Bird Report

hs.pop        <- seq(8, 80, 12)          # assume 4 eggs per clutch
                                        # i.e. 2 clutches to 20
                                        # lots of freedom to add to this
                                        # this is also going to vary massively 
                                        # based on spatial scale

capt.fail     <- 0.2                    # The number of HS birds that
                                        # fail during rearing & can't be
                                        # released

intro.years   <- seq(1, 5, 1)           # years in which hs birds will be
                                        # introduced

#----- Set phil.p and return.p to single values -----
# phil.p = l, m , h = 0.3, 0.5, 0.7
# return.p = l, h = 0.7, 0.9

phil.p        <- 0.3                    # probability of first years return
                                        # to release site

return.p      <- 0.9                    # probability of breeding birds 
                                        # returning to release site


# Potentially need to set the number of hs.pop as a specific value for each
#   simulation. Otherwise, still too many param combos.


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
            "data/NI_HS_scenario_A_parameters_lowphil_highfid.csv",
            sep = ",",
            row.names = FALSE)

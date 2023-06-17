# Loop through Extinction Cutoff Simulations for Plots

# Need to loop through the simulations to extract those for each pop
#   size.
#   Then loop through the extinction cut offs to obtain frequency 
#     of sims falling below these values

# Store results in long table
# years initN mean  lcl ucl

# set pop sizes
pop.size <- c("100", "150", "200", "250")

# set ext thresholds
ext.cutoff <- c(4, 20, 40, 80)

# create data frame for storage
NI.sim.out.df.long <- data.frame()
NI.sim.out.match.long <- data.frame()

sim.out.pop <- list()
sim.out.match <- list()

# loop
for (p in pop.size) {
  
  #extract simulations
  sim.out.split <- sim.out[names(NI.sim.out) == p]
  
  # then bind columnwise
  sim.out.split <- do.call(cbind, sim.out.split)
  
  # get means of rows
  sim.out.mean <- rowMeans(sim.out.split)
  
  # get ucl and lcl
  sim.out.split.lcl <- apply(sim.out.split, 1, quantile, probs = 0.1)
  sim.out.split.ucl <- apply(sim.out.split, 1, quantile, probs = 0.9)
  
  sim.out.pop[[p]] <- data.frame(mean  = sim.out.mean,
                                 lcl   = sim.out.split.lcl,
                                 ucl   = sim.out.split.ucl)
  
  for (t in ext.cutoff) {
    
    ext.match <- which(sim.out.split[31, ] < t,
                            arr.ind = TRUE)
    sim.out.match       <- sim.out.split[, ext.match]
    sim.out.match.mean  <- rowMeans(sim.out.match)
    sim.out.match.lcl   <- apply(sim.out.match, 1, quantile, probs = 0.1)
    sim.out.match.ucl   <- apply(sim.out.match, 1, quantile, probs = 0.9)
    
    sim.out.match[[t]]  <- data.frame(mean = sim.out.match.mean,
                                      lcl   = sim.out.match.lcl,
                                      ucl   = sim.out.match.ucl)
    
  }
  
  
}




bettermatch.100 <- which(NI.sim.out.100[31 ,] < 5,
                         arr.ind = TRUE)
NI.sim.out.100.match <- NI.sim.out.100[, bettermatch.100]
NI.sim.out.100.match.mean <- rowMeans(NI.sim.out.100.match)
NI.sim.out.100.match.lcl <- apply(NI.sim.out.100.match, 1,
                                  quantile, probs = 0.1)
NI.sim.out.100.match.ucl <- apply(NI.sim.out.100.match, 1,
                                  quantile, probs = 0.9)


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

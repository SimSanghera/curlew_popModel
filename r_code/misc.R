roi.5000 <- ROI.sim.out[names(ROI.sim.out) == "5000"]
roi.5000 <- do.call(cbind, roi.5000)
years <- seq(1:31)
row.names(roi.5000) <- years
roi.5000.mean <- rowMeans(roi.5000)
roi.lcl <- apply(roi.5000,1, quantile, probs=0.1)
roi.ucl <- apply(roi.5000, 1, quantile, probs = 0.9)

roi.5000.df <- data.frame(roi.5000.mean, roi.lcl, roi.ucl)
roi.5000.df$year <- seq(1:31)
colnames(roi.5000.df) <- c("pop.size", "lcl", "ucl", "year")

roi.5000.plot <- ggplot(roi.5000.df) + geom_line(aes(x = year,
                                                     y = pop.size)) + 
  geom_ribbon(aes(x = year,
                  ymin = lcl,
                  ymax = ucl),
              alpha = 0.3)

# set element names in list based on init.N
# extract and combine all columns with the same init.N
ROi.5000.df <- data.frame(years,
                          mean = ROI.sim.out.5000.mean,
                          lcl = ROI.sim.out.5000.lcl,
                          ucl = ROI.sim.out.5000.ucl)
str(ROi.5000.df)
roi.5000.plot <- ggplot(ROi.5000.df) +
  geom_line(aes(x = years,
                y = mean))
roi.5000.plot

roi.5000.ribbon <- roi.5000.plot + 
  geom_ribbon(aes(x = years,
                  ymin = lcl,
                  ymax = ucl),
              alpha = 0.3) +
  ylim(0, 15000)


vec1 <- NI.sim.out[[1]][, 1]
vec2 <- NI.sim.out[[1]][, 2]
vec3 <- NI.sim.out[[1]][, 3]
vec4 <- NI.sim.out[[2]][, 1]
vec5 <- NI.sim.out[[2]][, 2]
vec6 <- NI.sim.out[[2]][, 3]
vec7 <- NI.sim.out[[3]][, 1]
vec8 <- NI.sim.out[[3]][, 2]
vec9 <- NI.sim.out[[3]][, 3]
vec10 <- NI.sim.out[[4]][, 1]
vec11 <- NI.sim.out[[4]][, 2]
vec12 <- NI.sim.out[[4]][, 3]



mat1 <- as.matrix(vec1, vec2, vec3)

test.list <- list()
test.list[[1]] <- cbind(vec1, vec2, vec3)
test.list[[2]] <- cbind(vec4, vec5, vec6)
test.list[[3]] <- cbind(vec7, vec8, vec9)
test.list[[4]] <- cbind(vec10, vec11, vec12)

names(test.list) <- c("group1", "group2", "group1", "group2")
list.names <- (unique(names(test.list)))

test.list.2 <- sapply(unique(names(test.list)),
                      function(x) unname(cbind(test.list[names(test.list) == x])))





#-----  Plotting  -----

ROI.sim.out.plot <- ROI.sim.out %>%
  group_by()








proj.pops <- projections$pop.sizes
proj.pops.2 <- projections$pop.sizes
proj.pops.3 <- projections$pop.sizes

curlew.pop <- data.frame(proj.pops, proj.pops.2, proj.pops.3)
curlew.pop$mean <- rowMeans(curlew.pop)
curlew.pop$lcl <- apply(curlew.pop[1:3], 1, quantile, probs = 0.1)
curlew.pop$ucl <- apply(curlew.pop[1:3], 1, quantile, probs = 0.9)

curlew.pop$years <- years

curlew.plot <- ggplot(curlew.pop) +
  geom_line(aes(x = years,
                y = mean)) + 
  geom_ribbon(aes(x = years,
                  ymin = lcl,
                  ymax = ucl), alpha = 0.1)

curlew.plot.tidy <-curlew.plot +
  xlab("Year") +
  ylab("Population Size (individuals)") +
  theme(panel.background = element_rect(fill = "white", colour = "black"))

curlew.plot.tidy


curlew.pop.3800 <- data.frame(proj.pops, proj.pops.2, proj.pops.3)
curlew.pop.3800$mean <- rowMeans(curlew.pop.3800)
curlew.pop.3800$lcl <- apply(curlew.pop.3800[1:3], 1, quantile, probs = 0.1)
curlew.pop.3800$ucl <- apply(curlew.pop.3800[1:3], 1, quantile, probs = 0.9)

curlew.pop.3800$years <- years

curlew.plot.3800 <- ggplot(curlew.pop.3800) +
  geom_line(aes(x = years,
                y = mean)) + 
  geom_ribbon(aes(x = years,
                  ymin = lcl,
                  ymax = ucl), alpha = 0.5)

curlew.plot.3800.tidy <-curlew.plot.3800 +
  xlab("Year") +
  ylab("Population Size (individuals)") +
  theme(panel.background = element_rect(fill = "white", colour = "black"))

curlew.plot.3800.tidy
#-----  Curlew Code Edits and Trials  -----



  n.reps = 1          
  n.years = 30         
  init.N = (560 + 130)       
  K = 5000                
  catastrophe.p = 0.001   
  catastrophe.i = 0.25   
  F.loop = simul.in[50000, 4]         
  S.jw.loop = simul.in[50000, 3]        
  S.ad.loop = simul.in[50000, 2]      
  breed.fail = 0.1    

y = 2
rep = 1

Pop.Array <- array(0,
                   dim = c((n.years + 1),
                           n.reps))


Pop.Array[1, rep]   <- init.N*2

F.Year <- ifelse(rbinom(1, 1, breed.fail) == 1,
                 F.loop * 0.5, 
                 F.loop)

bird.vr <- list(prod.rates = F.Year,
                S.ad = S.ad.loop,
                S.jw = S.jw.loop #,
                # sx.r = sex.ratio
)

a <- matrix(sapply(bird.matrix,
                   eval,
                   bird.vr,
                   NULL),
            nrow = sqrt(length(bird.matrix)),
            byrow = TRUE)

pop.size  <- c(0, (init.N*2))

projections <- pop.projection(a,
                              n = pop.size,
                              iterations = 30)

R.max <- projections$lambda  

next.year <- max(0, 
                 trunc(Ricker(Pop.Array[y-1, rep])))

Pop.Array[y, rep] <- next.year







proj.pops <- projections$pop.sizes


years <- seq(2021, 2050, 1)
curlew.chicks <- as.vector(projections$stage.vectors[1 ,])
curlew.adults <- as.vector(projections$stage.vectors[2 ,])


curlew.pop <- as.data.frame(cbind(years, proj.pops))
curlew.pop$lambda <- c(projections$pop.changes, projections$lambda)
curlew.pop$chicks <- curlew.chicks
curlew.pop$adults <- curlew.adults

names(curlew.pop) <- c("year",
                       "total_pop",
                       "lambda",
                       "chicks",
                       "adults")

curlew.2050 <- ggplot(curlew.pop,
                      aes(x = years,
                          y = total_pop))

curlew.2050 <- curlew.2050 +
  geom_point() +
  geom_smooth()

curlew.2050 <- curlew.2050 +
  ylim(0, 1500) +
  xlab("Year") +
  ylab("Projected population size (individuals)") +
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

curlew.2050

curlew.splits.2050 <- curlew.2050 +
  geom_line(aes(y = chicks,
                x = years,
                colour = "darkred")) +
  geom_line(aes(y = adults,
                x = years,
                colour = "steelblue"))

curlew.splits.2050



a <- matrix(c(0, 0.2, 0, 0, 0,
              0, 0, 0.2, 0, 0,
              0.01, 0, 0, 0.326, 0,
              0.01, 0, 0, 0, 0.898,
              3, 0, 0, 0, 0.90),
            nrow = 5,
            ncol = 5,
            dimnames = list(c("egg",
                              "chick",
                              "First Yr",
                              "Second Yr",
                              "Adult"),
                            c("egg",
                              "chick",
                              "First Yr",
                              "Second Yr",
                              "Adult")))
a

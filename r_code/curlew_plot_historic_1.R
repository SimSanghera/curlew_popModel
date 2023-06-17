##### Curlew Plots  #####

#----- Cleaned script for plotting curlew model output  -----
#
# Libraries tidyverse, ggplot2
#
# Historic Model Plots
#   - The historic models were used to validate the model predictions.
#   - This fact is sometimes lost on the rest of the curlew brain trust.
#   - So you might need to keep hitting this on the head.
#   - But essentially, we ran the model - with fewer age classes - using 
#   the starting population from the 80s, then predicted trends to see
#   if they incorporated the "counts" from the 90s and 2000s and how closely
#   it predicted to these points.
#   - Counts is a massive misnomer, by the way. The faith in the curlew 
#   estimates from the 80s and 90s surveys is low. Especially for ROI.
#   - Initial population numbers were obtained from the literature, and using
#   the upper and lower bounds to cover the breadth of variation. 
#   - Whether this is important or not overall is debatable.
#   - ANYWAY, this script is just for the plots :-)
#   - So I decided to plot the different predictions in a grid form 
#   (essentially facet_wrap for initial pop)
#   - Adding points and error bars for each of the survey years
#   - Then an overall trend line
#
#--------------------------------------------------------------#


#----- Libraries  -----
libs <- c("tidyverse",
          "ggplot2",
          "fishualize",
          "ggpubr")

lapply(libs,
       require,
       character.only = TRUE)


#-----  Plots -----
NI.historic.df <- read.csv("data/NI_historic_1987_match750_1.csv",
                           sep = ",",
                           header = TRUE)

# Define inital starting populations & upper and lower limits
# Store in a data frame
pop <- c((5000), (2091), (526))
pop.lcl <- c((3800), (1243), (252))
pop.ucl <- c((6250), (2928), (783))
NI.pop.ests <- data.frame(years = c(1987, 1999, 2013),
                          popsize = pop,
                          pop.min = pop.lcl,
                          pop.max = pop.ucl)


NI.historic.plot <- ggplot(data = NI.historic.df) +
  # Set facet
  facet_wrap(~initN) +
  
  geom_line(aes(x = years,
                y = mean/2),
            size = 0.75) +
  
  geom_ribbon(aes(x = years,
                  ymin = lcl/2,
                  ymax = ucl/2),
              alpha = 0.1) +
  
  geom_point(data = NI.pop.ests,
             aes(x = years,
                 y = popsize),
             colour = "red") +
  geom_errorbar(data = NI.pop.ests,
                aes(x = years,
                    ymin = pop.min, 
                    ymax = pop.max),
                colour = "red",
                width = 1,
                size = 0.5) +
  
  # Set colour
  scale_colour_manual(values = megaPal2) +
  scale_fill_manual(values = megaPal2) +
  geom_vline(xintercept = 2034,
             linetype = "dashed", colour = "grey") +
  
  scale_y_continuous(
    expand = expansion(mult = c(0)),
    limits = c(0, 8000),
    breaks = seq(0, 8000, 500)) +
  scale_x_continuous(
    #  expand = c(0, 0),
    limits = c(1987, 2017),
    breaks = seq(1987, 2017, 3)
  ) +
  
  theme_bw() +
  theme(legend.position = c(0.7, 0.92),
        legend.title = element_text(size = 8,
                                    colour = "black",
                                    family = "Arial",
                                    face = "bold"),
        legend.text = element_text(size = 7),
        legend.direction = "horizontal") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(panel.grid = element_blank()) +
  
  theme(plot.title = element_blank()) +
  
  xlab("Years") +
  ylab("Population Size (breeding pairs)") +
  theme(axis.title.x = element_blank()) +
  theme(axis.ticks = element_line(size = 0.7)) +
  theme(axis.text.x = element_text(vjust = -0.5,
                                   size = 7,
                                   family = "Arial",
                                   angle = 0)) +
  theme(axis.text.y = element_text(size = 7,
                                   family = "Arial",
                                   margin = margin(r = 5))) +
  theme(axis.title.y = element_text(size = 10,
                                    family = "Arial",
                                    face = "bold",
                                    margin = margin(r = 5))) +
  theme(axis.title.x = element_text(size = 8,
                                    family = "Arial",
                                    face = "bold",
                                    margin = margin(t = 1,
                                                    r = 5,
                                                    b = 2),
                                    vjust = -1))

# Write to file
ggsave("figures/NI_historic_1.tiff",
       width = 9,
       height = 6,
       bg = "white") 


# RoI 
ROI.historic.df <- read.csv("data/ROI_sims_1985_match400.csv",
                            sep = ",",
                            header = TRUE)

# Define inital starting populations & upper and lower limits
# Store in a data frame
ROIpop <- c((7000), (1000), (138))
ROI.lcl <- c((5000), (743), (100))
ROI.ucl <- c((8500), (1850), (250))
ROI.pop.ests <- data.frame(years = c(1985, 2001, 2015),
                           popsize = ROIpop,
                           pop.min = ROI.lcl,
                           pop.max = ROI.ucl)

ROI.historic.plot <- ggplot(data = ROI.historic.df) +
  # Set facet
  facet_wrap(~initN) +
  
  geom_line(aes(x = years,
                y = mean/2),
            size = 0.75) +
  
  geom_ribbon(aes(x = years,
                  ymin = lcl/2,
                  ymax = ucl/2),
              alpha = 0.1) +
  
  geom_point(data = ROI.pop.ests,
             aes(x = years,
                 y = popsize),
             colour = "red") +
  geom_errorbar(data = ROI.pop.ests,
                aes(x = years,
                    ymin = pop.min, 
                    ymax = pop.max),
                colour = "red",
                width = 1,
                size = 0.5) +
  
  # Set colours
  scale_colour_manual(values = megaPal2) +
  scale_fill_manual(values = megaPal2) +
  geom_vline(xintercept = 2034,
             linetype = "dashed", colour = "grey") +
  
  scale_y_continuous(
    expand = expansion(mult = c(0)),
    limits = c(0, 9000),
    breaks = seq(0, 9000, 500)) +
  scale_x_continuous(
    #  expand = c(0, 0),
    limits = c(1985, 2015),
    breaks = seq(1985, 2015, 3)
  ) +
  
  theme_bw() +
  theme(legend.position = c(0.7, 0.92),
        legend.title = element_text(size = 8,
                                    colour = "black",
                                    family = "Arial",
                                    face = "bold"),
        legend.text = element_text(size = 7),
        legend.direction = "horizontal") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(panel.grid = element_blank()) +
  
  theme(plot.title = element_blank()) +
  
  xlab("Years") +
  ylab("Population Size (breeding pairs)") +
  theme(axis.title.x = element_blank()) +
  theme(axis.ticks = element_line(size = 0.7)) +
  theme(axis.text.x = element_text(vjust = -0.5,
                                   size = 7,
                                   family = "Arial",
                                   angle = 0)) +
  theme(axis.text.y = element_text(size = 7,
                                   family = "Arial",
                                   margin = margin(r = 5))) +
  theme(axis.title.y = element_text(size = 10,
                                    family = "Arial",
                                    face = "bold",
                                    margin = margin(r = 5))) +
  theme(axis.title.x = element_text(size = 8,
                                    family = "Arial",
                                    face = "bold",
                                    margin = margin(t = 1,
                                                    r = 5,
                                                    b = 2),
                                    vjust = -1))

# Write to file
ggsave("figures/ROI_historic_1.tiff",
       width = 9,
       height = 6,
       bg = "white") 



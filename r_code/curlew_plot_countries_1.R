##### Curlew Plots  #####

#----- Cleaned script for plotting curlew model output  -----
#
# Libraries tidyverse, ggplot2, fishualize, ggpubr
# Fishualize used for colour palettes. This can be replaced with base or 
#   any other colour palette of choice.
# ggpubr:: ggarrange() function for arranging the plots
#
# Country Level Plots:
#   - used the ribbon element to add upper and lower 
#   bands to the plot so we could visualise the uncertainty within the models,
#   with the mean plotted as a solid line.
# 
# The original model was created using 1000 simulations of far too many
#   parameters. These parameters were reduced following discussions with
#   Gillian, Anne-Marie and Steffen Oppel.
#   Given the range in values for many parameters, this still left a lot of
#   variation in model output - which we kept in deliberately to start off with.
#   If a full model is run and then plotted with ribbons, you will see the 
#   level of variation is high. 
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


#----- Define Colour Palettes -----
megaPal  <- fish(n = 8, option = "Lepomis_megalotis")
megaPal2 <- fish(n = 4, option = "Lepmois_megalotis")


#-----  NI Plots  -----

# Scenario A: 
#   - no change in demographic rates over time
#   - 
# NI has 3 plots for each scenario to cover the uncertainty of 
#   initial starting population size (150, 200, 250)


  #----- NI 150 -----
NI_A.df <- read.csv("data/curlew_HS_NI_A_max_simple_4.csv",
                    sep = ",",
                    header = T)

NI_A.plot <- ggplot(data = NI_A.df) +
  
  geom_line(aes(x      = years,
                y      = pop_size,
                group  = scenario,
                colour = scenario
                ),
            size = 0.75) +
  
  geom_ribbon(aes(x     = years,
                  ymin  = lcl,
                  ymax  = ucl,
                  group = scenario,
                  fill  = scenario),
              alpha = 0.1) +
  
  scale_colour_manual(values = megaPal2) +
  scale_fill_manual(values = megaPal2) +
  
  geom_vline(xintercept = 2034,
             linetype = "dashed",
             colour   = "grey") +
  
  scale_y_continuous(
    expand = expansion(mult = c(0)),
    limits = c(0, 400),
    breaks = seq(0, 400, 25)) + 
  
  scale_x_continuous(
    limits = c(2023, 2054),
    breaks = seq(2023, 2053, 1)) + 
  
  theme_bw() +
  theme(legend.position = c(0.7, 0.9),
        legend.title    = element_text(size   = 8,
                                       colour = "black",
                                       family = "Arial", # maybe remove family
                                       face   = "bold"),
        legend.text     = element_text(size   = 7),
        legend.direction = "horizontal") +
  
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(panel.grid       = element_blank(),
        panel.grid       = element_blank()) +
  
  theme(plot.title = element_blank()) +
  
  xlab("Years") +
  ylab("Population Size (breeding pairs)") +
  theme(axis.title.x = element_blank()) +
  theme(axis.ticks   = element_line(size   = 0.7)) +
  theme(axis.text.x  = element_text(vjust  = -0.5,
                                    size   = 7,
                                    family = "Arial",
                                    angle  = 0)) +
  theme(axis.text.y  = element_text(size   = 7,
                                    family = "Arial",
                                    margin = margin(r = 5))) +
  theme(axis.title.y = element_text(size   = 10,
                                    family = "Arial",
                                    face   = "bold",
                                    margin = margin(t = 1,
                                                    r = 5,
                                                    b = 2),
                                    vjust  = -1))

# Write to file
ggsave("figures/NI_A_max_150_ci_2.tiff",
       width = 9,
       height = 6,
       bg = "white")


  #----- NI 200 -----
NI_A_200.df <- read.csv("data/curlew_HS_NI_A200_max_simple_1.csv",
                        sep = ",",
                        header = T)
NI_A_200.plot <- ggplot(data = NI_A_200.df) +
  
  geom_line(aes(x      = years,
                y      = pop_size,
                group  = scenario,
                colour = scenario
                ),
            size = 0.75) +
  
  geom_ribbon(aes(x     = years,
                  ymin  = lcl,
                  ymax  = ucl,
                  group = scenario,
                  fill  = scenario),
              alpha = 0.1) +
  
  scale_colour_manual(values = megaPal2) +
  scale_fill_manual(values = megaPal2) +
  
  geom_vline(xintercept = 2034,
             linetype = "dashed",
             colour   = "grey") +
  
  theme_bw() +
  theme(legend.position = c(0.7, 0.9),
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
  theme(axis.title.x = element_text(size = 12,
                                    family = "Arial",
                                    face = "bold",
                                    margin = margin(t = 1,
                                                    r = 5,
                                                    b = 2),
                                    vjust = -1))

# Write to file
ggsave("figures/NI_A_max_200_1.tiff",
       width = 9,
       height = 6,
       bg = "white") 



  #----- NI 250 -----
NI_A_250.df <- read.csv("data/curlew_HS_NI_A250_max_simple_1.csv",
                        sep = ",",
                        header = T)

NI_A_250.plot <- ggplot(data = NI_A_250.df) +
  
  geom_line(aes(x      = years,
                y      = pop_size,
                group  = scenario,
                colour = scenario),
            size = 0.75) +
  
  geom_ribbon(aes(x     = years,
                  ymin  = lcl,
                  ymax  = ucl,
                  group = scenario,
                  fill  = scenario),
              alpha = 0.1) +
  
  scale_colour_manual(values = megaPal2) +
  scale_fill_manual(values = megaPal2) +
  
  geom_vline(xintercept = 2034,
             linetype = "dashed",
             colour   = "grey") +
  
  scale_colour_manual(values = megaPal2) +
  scale_fill_manual(values = megaPal2) +
  geom_vline(xintercept = 2034,
             linetype = "dashed", colour = "grey") +
  
  scale_y_continuous(
    expand = expansion(mult = c(0)),
    limits = c(0, 500),
    breaks = seq(0, 500, 25)) +
  scale_x_continuous(
    #  expand = c(0, 0),
    limits = c(2023, 2054),
    breaks = seq(2023, 2053, 1)
  ) +
  
  theme_bw() +
  theme(legend.position = c(0.7, 0.9),
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
  theme(axis.title.x = element_text(size = 12,
                                    family = "Arial",
                                    face = "bold",
                                    margin = margin(t = 1,
                                                    r = 5,
                                                    b = 2),
                                    vjust = -1))

# Write to file
ggsave("figures/NI_A_max_250_ci_1.tiff",
       width = 9,
       height = 6,
       bg = "white") 

  

#----- NI Scenario B  -----
NI_B.df <- read.csv("data/HS_NI_B_means_simple_2.csv",
                    sep = ",",
                    header = T)

NI_B.plot <- ggplot(data = NI_B.df) +
  
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            size = 0.75) +
  geom_ribbon(aes(x = years,
                  ymin = lcl,
                  ymax = ucl,
                  group = scenario,
                  fill = scenario),
              alpha = 0.1) +
  
  
  scale_colour_manual(values = megaPal2) +
  scale_fill_manual(values = megaPal2) +
  geom_vline(xintercept = 2034,
             linetype = "dashed", colour = "grey") +
  
  scale_y_continuous(
    expand = expansion(mult = c(0)),
    limits = c(0, 400),
    breaks = seq(0, 400, 25)) +
  scale_x_continuous(
    #  expand = c(0, 0),
    limits = c(2023, 2054),
    breaks = seq(2023, 2053, 1)
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
  theme(axis.title.x = element_text(size = 12,
                                    family = "Arial",
                                    face = "bold",
                                    margin = margin(t = 1,
                                                    r = 5,
                                                    b = 2),
                                    vjust = -1))

# Write to file
ggsave("figures/NI_B_max_150_6.tiff",
       width = 9,
       height = 6,
       bg = "white") 





#-----  ROI -----
# Scenario A: 
#   - no change in demographic rates over time
#   - 

ROI_A.df <- read.csv("data/HS_ROI_A_means_simple_3.csv",
                     sep = ",",
                     header = T)

ROI_A.plot <- ggplot(data = ROI_A.df) +

    geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            size = 0.75) +
  geom_ribbon(aes(x = years,
                  ymin = lcl,
                  ymax = ucl,
                  group = scenario,
                  fill = scenario),
              alpha = 0.1) +

  scale_colour_manual(values = megaPal2) +
  scale_fill_manual(values = megaPal2) +
  geom_vline(xintercept = 2034,
             linetype = "dashed", colour = "grey") +
  
  scale_y_continuous(
    expand = expansion(mult = c(0)),
    limits = c(0, 200),
    breaks = seq(0, 200, 25)) +
  scale_x_continuous(
    #  expand = c(0, 0),
    limits = c(2023, 2054),
    breaks = seq(2023, 2053, 1)
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
  theme(axis.title.x = element_text(size = 12,
                                    family = "Arial",
                                    face = "bold",
                                    margin = margin(t = 1,
                                                    r = 5,
                                                    b = 2),
                                    vjust = -1))

# Write to file
ggsave("figures/ROI_A_max_6.tiff",
       width = 9,
       height = 6,
       bg = "white") 


#----- Scenario B -----
ROI_B.df <- read.csv("data/HS_ROI_B_simple_3.csv",
                     sep = ",",
                     header = T)

ROI_B.plot <- ggplot(data = ROI_B.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            size = 0.75) +
  geom_ribbon(aes(x = years,
                  ymin = lcl,
                  ymax = ucl,
                  group = scenario,
                  fill = scenario),
              alpha = 0.1) +
  
  
  scale_colour_manual(values = megaPal2) +
  scale_fill_manual(values = megaPal2) +
  geom_vline(xintercept = 2034,
             linetype = "dashed", colour = "grey") +
  
  scale_y_continuous(
    expand = expansion(mult = c(0)),
    limits = c(0, 250),
    breaks = seq(0, 250, 25)) +
  scale_x_continuous(
    #  expand = c(0, 0),
    limits = c(2023, 2054),
    breaks = seq(2023, 2053, 1)
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
  theme(axis.title.x = element_text(size = 12,
                                    family = "Arial",
                                    face = "bold",
                                    margin = margin(t = 1,
                                                    r = 5,
                                                    b = 2),
                                    vjust = -1))

# Write to file
ggsave("figures/ROI_B_max_ci_2.tiff",
       width = 9,
       height = 6,
       bg = "white") 




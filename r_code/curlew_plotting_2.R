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
# Site Level Plots:
#   - for expediency and ease-of-reading for non-modellers, the
#   site level plots were simply plotted using single lines for each scenario
#   that was run (number of eggs introduced into the system)
#   - the plots became overly cluttered when introducing upper and lower CIs
#
#--------------------------------------------------------------#



#----- Load Libraries -----
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




#----- Plots for Elasticity -----
# This was simply to visualise the impact of each demographic paramater on
#   the predicted lambda (growth rate).
elas.df <- read.csv("data/elasticity_analysis.csv",
                    sep = ",",
                    header = T)

elas.plot <- ggplot(data = elas.df,
                    aes(x = Parameters,
                        y = Lambda,
                        group = Values,
                        fill  = Values)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  geom_text(aes(group = Values,
                label = round(PercDiff, digits = 2)),
            position = position_dodge(.9),
            size = 2,
            vjust = -1) +
  scale_fill_manual(values = kisPal) +
  geom_hline(yintercept = 0.7457557,
             linetype = "dashed", colour = "grey") +
  
  scale_y_continuous(expand = expansion(mult = c(0)),
                     limits = c(0, 1),
                     breaks = seq(0, 1, 0.2)) +
  
  theme_bw() +
  theme(legend.position = c(0.95, 0.91),
        legend.title = element_text(size = 7,
                                    colour = "black",
                                    family = "Arial"),
        legend.text = element_text(size = 7)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(panel.grid = element_blank()) +
  
  theme(plot.title = element_blank()) +
  
  ylab("Population Growth Rate") +
  xlab("Parameters") +
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

# Write plot to file
ggsave("figures/elasticity_analysis.tiff",
       width = 9,
       height = 6,
       bg = "white")  

#---------------------------------------------------#



#----- NI Site Plots -----
# Glenwherry & Lower Loch Erne (LLE)

  #---- Glenwherry  -----
glen.df <- read.csv("data/curlew_glen_20bp_means_simple_2.csv",
                    sep = ",",
                    header = T)

# change labels in data
glen.df$num_clutches <- ifelse(glen.df$num_clutches == "base",
                               "current_state",
                               glen.df$num_clutches)

glen.plot <- ggplot(data = glen.df,
                    aes(x = years,
                        y = pop_size,
                        group  = num_clutches,
                        colour = num_clutches),
                    size = 0.75) + 
  
  scale_colour_manual(values = megaPal) + # set colour palette
  geom_vline(xintercept = 2034,           # add vertical line for 10 year step
             linetype = "dashed", colour = "grey") +
  
  scale_y_continuous(                     # alter scale
    expand = expansion(mult = c(0)),
    limits = c(0, 110),
    breaks = seq(0, 110, 10)) +
  scale_x_continuous(
    #  expand = c(0, 0),
    limits = c(2023, 2054),
    breaks = seq(2023, 2053, 1)
  ) +
  
  # Maybe required to remove "family" paramater from code as it doesn't seem
  #   to have any effect
  theme_bw() +
  theme(legend.position = c(0.7, 0.15),
        legend.title = element_text(size = 8,
                                    colour = "black",
                                    family = "Arial",
                                    face = "bold"),
        legend.text = element_text(size = 7),
        legend.direction = ("horizontal")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(panel.grid = element_blank()) +
  theme(panel.grid = element_blank()) +
  
  theme(plot.title = element_blank()) +
  
  xlab("Years") +
  ylab("Population Size (breeding pairs)") +
  labs(colour = expression("Scenario")) +
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


# Write to File
ggsave("figures/glenwherry_4.tiff",
       width = 9,
       height = 6,
       bg = "white") 


  #----- LLE  -----
lle.df <- read.csv("data/curlew_lle_20bp_means_simple_2.csv",
                   sep = ",",
                   header = T)
# change labels
lle.df$num_clutches <- ifelse(lle.df$num_clutches == "base",
                              "current_state",
                              lle.df$num_clutches)

lle.plot <- ggplot(data = lle.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = num_clutches,
                colour = num_clutches),
            size = 0.75) +

  scale_colour_manual(values = megaPal) +
  geom_vline(xintercept = 2034,
             linetype = "dashed", colour = "grey") +
  
  scale_y_continuous(
    expand = expansion(mult = c(0)),
    limits = c(0, 90),
    breaks = seq(0, 90, 10)) +
  scale_x_continuous(
    #  expand = c(0, 0),
    limits = c(2023, 2054),
    breaks = seq(2023, 2053, 1)
  ) +
  
  theme_bw() +
  theme(legend.position = c(0.7, 0.15),
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
  labs(colour = expression("Scenario")) +
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
ggsave("figures/lle_3.tiff",
       width = 9,
       height = 6,
       bg = "white") 


#--------------------------------------------------#






#----- ROI Sites  -----
# Corrib, Stacks, Donegal, Monaghan, Ree  

  #---- Corrib  -----
corrib.df <- read.csv("data/curlew_HS_corrib_means_simple_2.csv",
                      sep = ",",
                      header = T)

# Change labels
corrib.df$scenario <- ifelse(corrib.df$scenario == "base", "current_state",
                             corrib.df$scenario)
corrib.df$scenario <- ifelse(corrib.df$scenario == "HS", "feasible_clutches",
                             corrib.df$scenario)
corrib.df$scenario <- ifelse(corrib.df$scenario == "HS_max", "max_clutches",
                             corrib.df$scenario)


corrib.plot <- ggplot(data = corrib.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            size = 0.75) +
  
  scale_colour_manual(values = megaPal2) +
  geom_vline(xintercept = 2034,
             linetype = "dashed", colour = "grey") +
  
  scale_y_continuous(
    expand = expansion(mult = c(0)),
    limits = c(0, 40),
    breaks = seq(0, 40, 10)) +
  scale_x_continuous(
    #  expand = c(0, 0),
    limits = c(2023, 2054),
    breaks = seq(2023, 2053, 1)
  ) +
  
  theme_bw() +
  theme(legend.position = c(0.7, 0.15),
        legend.title = element_text(size = 8,
                                    colour = "black",
                                    family = "Arial",
                                    face = "bold"),
        legend.text = element_text(size = 7),
        legend.direction = "horizontal") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
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

# Write to File
ggsave("figures/corrib_3.tiff",
       width = 9,
       height = 6,
       bg = "white") 


  #-----  Stacks  -----
stacks.df <- read.csv("data/curlew_HS_stacks_means_simple_1.csv",
                      sep = ",",
                      header = T)
# Change labels
stacks.df$scenario <- ifelse(stacks.df$scenario == "base", "current_state",
                             stacks.df$scenario)
stacks.df$scenario <- ifelse(stacks.df$scenario == "HS", "feasible_clutches",
                             stacks.df$scenario)
stacks.df$scenario <- ifelse(stacks.df$scenario == "HS_max", "max_clutches",
                             stacks.df$scenario)

stacks.plot <- ggplot(data = stacks.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            size = 0.75) +
  
  scale_colour_manual(values = megaPal2) +
  geom_vline(xintercept = 2034,
             linetype = "dashed", colour = "grey") +
  
  scale_y_continuous(
    expand = expansion(mult = c(0)),
    limits = c(0, 10),
    breaks = seq(0, 10, 1)) +
  scale_x_continuous(
    #  expand = c(0, 0),
    limits = c(2023, 2054),
    breaks = seq(2023, 2053, 1)
  ) +
  
  theme_bw() +
  theme(legend.position = c(0.7, 0.35),
        legend.title = element_text(size = 7,
                                    colour = "black",
                                    family = "Arial"),
        legend.text = element_text(size = 7),
        legend.direction = "horizontal") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
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
ggsave("figures/stacks_3.tiff",
       width = 9,
       height = 6,
       bg = "white") 


  #----- Ree  -----
ree.df <- read.csv("data/curlew_HS_ree_means_simple_2.csv",
                   sep = ",",
                   header = T)
# Change labels
ree.df$scenario <- ifelse(ree.df$scenario == "base", "current_state",
                          corrib.df$scenario)
ree.df$scenario <- ifelse(ree.df$scenario == "HS", "feasible_clutches",
                          ree.df$scenario)
ree.df$scenario <- ifelse(ree.df$scenario == "HS_max", "max_clutches",
                          ree.df$scenario)

ree.plot <- ggplot(data = ree.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            position = position_dodge(width = 0.5),
            size = 0.75) +
  
  scale_colour_manual(values = megaPal2) +
  geom_vline(xintercept = 2034,
             linetype = "dashed", colour = "grey") +
  
  scale_y_continuous(
    expand = expansion(mult = c(0)),
    limits = c(0, 30),
    breaks = seq(0, 30, 5)) +
  scale_x_continuous(
    #  expand = c(0, 0),
    limits = c(2023, 2054),
    breaks = seq(2023, 2053, 1)
  ) +
  
  theme_bw() +
  theme(legend.position = c(0.7, 0.15),
        legend.title = element_text(size = 8,
                                    colour = "black",
                                    family = "Arial",
                                    face = "bold"),
        legend.text = element_text(size = 7),
        legend.direction = "horizontal") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
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
ggsave("figures/ree_2.tiff",
       width = 9,
       height = 6,
       bg = "white") 


  #----- Monaghan -----
mon.df <- read.csv("data/curlew_HS_mon_means_simple_1.csv",
                   sep = ",",
                   header = T)
# Change labels
mon.df$scenario <- ifelse(mon.df$scenario == "base", "current_state",
                          mon.df$scenario)
mon.df$scenario <- ifelse(mon.df$scenario == "HS", "feasible_clutches",
                          mon.df$scenario)
mon.df$scenario <- ifelse(mon.df$scenario == "HS_max", "max_clutches",
                          mon.df$scenario)

mon.plot <- ggplot(data = mon.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            size = 0.75) +
  
  scale_colour_manual(values = megaPal2) +
  geom_vline(xintercept = 2034,
             linetype = "dashed", colour = "grey") +
  
  scale_y_continuous(
    expand = expansion(mult = c(0)),
    limits = c(0, 10),
    breaks = seq(0, 10, 1)) +
  scale_x_continuous(
    #  expand = c(0, 0),
    limits = c(2023, 2054),
    breaks = seq(2023, 2053, 1)
  ) +
  
  theme_bw() +
  theme(legend.position = c(0.7, 0.82),
        legend.title = element_text(size = 7,
                                    colour = "black",
                                    family = "Arial"),
        legend.text = element_text(size = 7),
        legend.direction = "horizontal") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
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
ggsave("figures/mon_2.tiff",
       width = 9,
       height = 6,
       bg = "white") 


  #----- Donegal -----
don.df <- read.csv("data/curlew_HS_don_means_simple_1.csv",
                   sep = ",",
                   header = T)
# Change labels
don.df$scenario <- ifelse(don.df$scenario == "base", "current_state",
                          don.df$scenario)
don.df$scenario <- ifelse(don.df$scenario == "HS", "feasible_clutches",
                          don.df$scenario)
don.df$scenario <- ifelse(don.df$scenario == "HS_max", "max_clutches",
                          don.df$scenario)

don.plot <- ggplot(data = don.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            size = 0.75) +
  
  scale_colour_manual(values = megaPal2) +
  geom_vline(xintercept = 2034,
             linetype = "dashed", colour = "grey") +
  
  scale_y_continuous(
    expand = expansion(mult = c(0)),
    limits = c(0, 10),
    breaks = seq(0, 10, 1)) +
  scale_x_continuous(
    #  expand = c(0, 0),
    limits = c(2023, 2054),
    breaks = seq(2023, 2053, 1)
  ) +
  
  theme_bw() +
  theme(legend.position = c(0.7, 0.82),
        legend.title = element_text(size = 7,
                                    colour = "black",
                                    family = "Arial"),
        legend.text = element_text(size = 7),
        legend.direction = "horizontal") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
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
ggsave("figures/don_2.tiff",
       width = 9,
       height = 6,
       bg = "white") 


#---------------------------------------------------------------#

#----- Plots of Removal Scenarios -----
#
# ROI: Removing set numbers of clutches (see report)
#   - 2, 4, 8, 12, 20
# First create the individual plots
# Then create the arrangement of all plots together


  #----- Removing 2 clutches  -----
rem2.df <- read.csv("data/curlew_HS_rem2_means_simple_1.csv",
                    sep = ",",
                    header = T)

rem2.plot <- ggplot(data = rem2.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            size = 0.75) +
  
  # Set colour palette - 4 colours
  scale_colour_manual(values = megaPal) +
  geom_vline(xintercept = 2034,
             linetype = "dashed", colour = "grey") +
  
  scale_y_continuous(
    expand = expansion(mult = c(0)),
    limits = c(0, 10),
    breaks = seq(0, 10, 1)) +
  scale_x_continuous(
    #  expand = c(0, 0),
    limits = c(2023, 2054),
    breaks = seq(2023, 2053, 1)
  ) +
  
  theme_bw() +
  theme(legend.position = c(0.7, 0.82),
        legend.title = element_text(size = 10,
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
ggsave("figures/rem2.tiff",
       width = 9,
       height = 6,
       bg = "white") 


  #----- Removing 4 clutches  -----
rem4.df <- read.csv("data/curlew_HS_rem4_means_simple_1.csv",
                    sep = ",",
                    header = T)

rem4.plot <- ggplot(data = rem4.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            size = 0.75) +
  
  scale_colour_manual(values = megaPal) +
  geom_vline(xintercept = 2034,
             linetype = "dashed", colour = "grey") +
  
  scale_y_continuous(
    expand = expansion(mult = c(0)),
    limits = c(0, 10),
    breaks = seq(0, 10, 1)) +
  scale_x_continuous(
    #  expand = c(0, 0),
    limits = c(2023, 2054),
    breaks = seq(2023, 2053, 1)
  ) +
  
  theme_bw() +
  theme(legend.position = c(0.7, 0.82),
        legend.title = element_text(size = 10,
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
ggsave("figures/rem4.tiff",
       width = 9,
       height = 6,
       bg = "white") 


  #----- Removing 8 Clutches  -----
rem8.df <- read.csv("data/curlew_HS_rem8_means_simple_1.csv",
                    sep = ",",
                    header = T)

rem8.plot <- ggplot(data = rem8.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            size = 0.75) +
  
  scale_colour_manual(values = megaPal) +
  geom_vline(xintercept = 2034,
             linetype = "dashed", colour = "grey") +
  
  scale_y_continuous(
    expand = expansion(mult = c(0)),
    limits = c(0, 10),
    breaks = seq(0, 10, 1)) +
  scale_x_continuous(
    #  expand = c(0, 0),
    limits = c(2023, 2054),
    breaks = seq(2023, 2053, 1)
  ) +
  
  theme_bw() +
  theme(legend.position = c(0.7, 0.82),
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
ggsave("figures/rem8.tiff",
       width = 9,
       height = 6,
       bg = "white") 


  #-----  Removing 12 clutches  -----
rem12.df <- read.csv("data/curlew_HS_rem12_means_simple_1.csv",
                     sep = ",",
                     header = T)

rem12.plot <- ggplot(data = rem12.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            size = 0.75) +
  
  scale_colour_manual(values = megaPal) +
  geom_vline(xintercept = 2034,
             linetype = "dashed", colour = "grey") +
  
  scale_y_continuous(
    expand = expansion(mult = c(0)),
    limits = c(0, 20),
    breaks = seq(0, 20, 2)) +
  scale_x_continuous(
    #  expand = c(0, 0),
    limits = c(2023, 2054),
    breaks = seq(2023, 2053, 1)
  ) +
  
  theme_bw() +
  theme(legend.position = c(0.7, 0.82),
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
ggsave("figures/rem12.tiff",
       width = 9,
       height = 6,
       bg = "white")


  #-----  Removing 20 Clutches  -----
rem20.df <- read.csv("data/curlew_HS_rem20_means_simple_1.csv",
                     sep = ",",
                     header = T)

rem20.plot <- ggplot(data = rem20.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            size = 0.75) +
  
  scale_colour_manual(values = megaPal) +
  geom_vline(xintercept = 2034,
             linetype = "dashed", colour = "grey") +
  
  scale_y_continuous(
    expand = expansion(mult = c(0)),
    limits = c(0, 30),
    breaks = seq(0, 30, 2)) +
  scale_x_continuous(
    #  expand = c(0, 0),
    limits = c(2023, 2054),
    breaks = seq(2023, 2053, 1)
  ) +
  
  theme_bw() +
  theme(legend.position = c(0.7, 0.82),
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
ggsave("figures/rem20.tiff",
       width = 9,
       height = 6,
       bg = "white") 



#----- Combine plots in grid  -----
# Add labels
# 2 by 3 grid

allrem <- ggarrange(
  rem2.plot,
  rem4.plot,
  rem8.plot,
  rem12.plot,
  rem20.plot,
  common.legend = TRUE,
  legend        = "bottom",
  labels        = c("A", "B", "C", "D", "E"),
  ncol          = 2, nrow = 3
)





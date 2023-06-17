##################

#----- Plots -----

fishualize(n = 8, option = "Epinephelus_fasciatus")
fishualize(n = 8, option = "Trimma_lantana")
fishualize(n = 4, option = "Oncorhynchus_kisutch")
fishualize(n = 8, option = "Etheostoma_spectabile")
fishualize(n = 8, option = "Cirrhilabrus_solorensis")
fishualize(n = 4, option = "Lepomis_megalotis")

fasPal <- fish(n = 8, option = "Epinephelus_fasciatus")
kisPal <- fish(n = 4, option = "Oncorhynchus_kisutch")
specPal <- fish(n = 8, option = "Etheostoma_spectabile")
soloPal <- fish(n = 8, option = "Cirrhilabrus_solorensis")
megaPal <- fish(n = 8, option = "Lepomis_megalotis")
megaPal2 <- fish(n = 4, option = "Lepomis_megalotis")


# Elasticity
elas.df <- read.csv("data/elasticity_analysis.csv",
                    sep = ",",
                    header = T)
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

ggsave("figures/elasticity_analysis.tiff",
       width = 9,
       height = 6,
       bg = "white")  



# scale_y_continuous(expand = expansion(mult = c(0)),
#                    limits = c(0.0, 1.0),
#                    breaks = seq(0.0, 1.0, 0.1)) +
# coord_cartesian(ylim = c(0.0, 1.0)) + 
# labs(colour = "Parameters") +
# theme_bw() +
# theme(axis.title = element_text(size = 32),
#       axis.title.x = element_text(vjust = -0.5),
#       axis.text = element_text(size = 16, colour = "black"),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       panel.border = element_rect(),
#       legend.title = element_text(size = 20),
#       legend.position = c(0.95, 0.91),
#       legend.justification = "right",
#       legend.text = element_text(size = 20),
#       strip.text.x = element_text(size = 20))


#----- Sites -----
# NI: glen, lle
# glen = curlew_glen_20bp_means_simple_2
glen.df <- read.csv("data/curlew_glen_20bp_means_simple_2.csv",
                    sep = ",",
                    header = T)
glen.df$num_clutches <- ifelse(glen.df$num_clutches == "base", "current_state",
                               glen.df$num_clutches)

glen.plot <- ggplot(data = glen.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = num_clutches,
                colour = num_clutches),
            size = 0.75) +
  # geom_smooth(aes(x = years,
  #                 y = pop_size,
  #                 group = num_clutches,
  #                 colour = num_clutches,
  #                 fill = num_clutches),
  #             size = 0.25,
  #             span = 0.3) +

  scale_colour_manual(values = megaPal) +
  geom_vline(xintercept = 2034,
             linetype = "dashed", colour = "grey") +

  scale_y_continuous(
    expand = expansion(mult = c(0)),
    limits = c(0, 110),
    breaks = seq(0, 110, 10)) +
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

ggsave("figures/glenwherry_4.tiff",
       width = 9,
       height = 6,
       bg = "white") 


# lle - curlew_lle_20bp_means_simple_2
lle.df <- read.csv("data/curlew_lle_20bp_means_simple_2.csv",
                    sep = ",",
                    header = T)
lle.df$num_clutches <- ifelse(lle.df$num_clutches == "base", "current_state",
                               lle.df$num_clutches)

lle.plot <- ggplot(data = lle.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = num_clutches,
                colour = num_clutches),
            size = 0.75) +
  # geom_smooth(aes(x = years,
  #                 y = pop_size,
  #                 group = num_clutches,
  #                 colour = num_clutches,
  #                 fill = num_clutches),
  #             size = 0.25,
  #             span = 0.3) +
  
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

ggsave("figures/lle_3.tiff",
       width = 9,
       height = 6,
       bg = "white") 

# RoI: stacks, corrib, don, mon, ree
# curlew_HS_corrib_means_simple_1
corrib.df <- read.csv("data/curlew_HS_corrib_means_simple_2.csv",
                   sep = ",",
                   header = T)
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
  # geom_smooth(aes(x = years,
  #                 y = pop_size,
  #                 group = scenario,
  #                 colour = scenario,
  #                 fill = scenario),
  #             size = 0.25,
  #             span = 0.3) +
  
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

ggsave("figures/corrib_3.tiff",
       width = 9,
       height = 6,
       bg = "white") 

# curlew_HS_stacks_means_simple_1
stacks.df <- read.csv("data/curlew_HS_stacks_means_simple_1.csv",
                      sep = ",",
                      header = T)
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
  # geom_smooth(aes(x = years,
  #                 y = pop_size,
  #                 group = scenario,
  #                 colour = scenario,
  #                 fill = scenario),
  #             size = 0.25,
  #             span = 0.3) +
  # 
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

ggsave("figures/stacks_3.tiff",
       width = 9,
       height = 6,
       bg = "white") 

# curlew_HS_ree_means_simple_1
ree.df <- read.csv("data/curlew_HS_ree_means_simple_2.csv",
                      sep = ",",
                      header = T)
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
  # geom_smooth(aes(x = years,
  #                 y = pop_size,
  #                 group = scenario,
  #                 colour = scenario,
  #                 fill = scenario),
  #             size = 0.25,
  #             span = 0.3) +
  
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

ggsave("figures/ree_2.tiff",
       width = 9,
       height = 6,
       bg = "white") 

# curlew_HS_mon_means_simple_1
mon.df <- read.csv("data/curlew_HS_mon_means_simple_1.csv",
                      sep = ",",
                      header = T)
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
  # geom_smooth(aes(x = years,
  #                 y = pop_size,
  #                 group = scenario,
  #                 colour = scenario,
  #                 fill = scenario),
  #             size = 0.25,
  #             span = 0.3) +
  
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

ggsave("figures/mon_2.tiff",
       width = 9,
       height = 6,
       bg = "white") 

# curlew_HS_don_means_simple_1
don.df <- read.csv("data/curlew_HS_don_means_simple_1.csv",
                      sep = ",",
                      header = T)
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
  # geom_smooth(aes(x = years,
  #                 y = pop_size,
  #                 group = scenario,
  #                 colour = scenario,
  #                 fill = scenario),
  #             size = 0.25,
  #             span = 0.3) +
  
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

ggsave("figures/don_2.tiff",
       width = 9,
       height = 6,
       bg = "white") 





#----- Removals -----

# RoI: 2, 4, 8, 12, 20
# curlew_HS_rem2_means_simple_1
rem2.df <- read.csv("data/curlew_HS_rem2_means_simple_1.csv",
                   sep = ",",
                   header = T)
rem2.plot <- ggplot(data = rem2.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            size = 0.75) +
  # geom_smooth(aes(x = years,
  #                 y = pop_size,
  #                 group = scenario,
  #                 colour = scenario,
  #                 fill = scenario),
  #             size = 0.25,
  #             span = 0.3) +
  
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

ggsave("figures/rem2.tiff",
       width = 9,
       height = 6,
       bg = "white") 

# curlew_HS_rem4_means_simple_1
rem4.df <- read.csv("data/curlew_HS_rem4_means_simple_1.csv",
                   sep = ",",
                   header = T)
rem4.plot <- ggplot(data = rem4.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            size = 0.75) +
  # geom_smooth(aes(x = years,
  #                 y = pop_size,
  #                 group = scenario,
  #                 colour = scenario,
  #                 fill = scenario),
  #             size = 0.25,
  #             span = 0.3) +
  
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

ggsave("figures/rem4.tiff",
       width = 9,
       height = 6,
       bg = "white") 

# curlew_HS_rem8_means_simple_1
rem8.df <- read.csv("data/curlew_HS_rem8_means_simple_1.csv",
                   sep = ",",
                   header = T)
rem8.plot <- ggplot(data = rem8.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            size = 0.75) +
  # geom_smooth(aes(x = years,
  #                 y = pop_size,
  #                 group = scenario,
  #                 colour = scenario,
  #                 fill = scenario),
  #             size = 0.25,
  #             span = 0.3) +
  
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

ggsave("figures/rem8.tiff",
       width = 9,
       height = 6,
       bg = "white") 

# curlew_HS_rem12_means_simple_1
rem12.df <- read.csv("data/curlew_HS_rem12_means_simple_1.csv",
                   sep = ",",
                   header = T)
rem12.plot <- ggplot(data = rem12.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            size = 0.75) +
  # geom_smooth(aes(x = years,
  #                 y = pop_size,
  #                 group = scenario,
  #                 colour = scenario,
  #                 fill = scenario),
  #             size = 0.25,
  #             span = 0.3) +
  
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

ggsave("figures/rem12.tiff",
       width = 9,
       height = 6,
       bg = "white") 

# curlew_HS_rem20_means_simple_1
rem20.df <- read.csv("data/curlew_HS_rem20_means_simple_1.csv",
                   sep = ",",
                   header = T)
rem20.plot <- ggplot(data = rem20.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            size = 0.75) +
  # geom_smooth(aes(x = years,
  #                 y = pop_size,
  #                 group = scenario,
  #                 colour = scenario,
  #                 fill = scenario),
  #             size = 0.25,
  #             span = 0.3) +
  
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

ggsave("figures/rem20.tiff",
       width = 9,
       height = 6,
       bg = "white") 


plot.list <- 
  

ggarrange(rem4.plot,
          rem8.plot,
          labels = c("A", "B"))    
  
allrem <- ggarrange(rem2.plot,
          rem4.plot,
          rem8.plot,
          rem12.plot,
          rem20.plot,
          common.legend = TRUE,
          legend = "bottom",
          labels = c("A", "B", "C", "D", "E"),
          ncol = 2, nrow = 3)

ggsave("figures/allrem.tiff",
       width = 18,
       height = 12,
       bg = "white")



#----- NI ------
# Scenario A all HS_NI_A_means_simple_1
NI_A.df <- read.csv("data/curlew_HS_NI_A_max_simple_4.csv",
                   sep = ",",
                   header = T)
NI_A.plot <- ggplot(data = NI_A.df) +
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

ggsave("figures/NI_A_max_150_ci_2.tiff",
       width = 9,
       height = 6,
       bg = "white") 

# Scenario A all NI 200 pop
NI_A_200.df <- read.csv("data/curlew_HS_NI_A200_max_simple_1.csv",
                    sep = ",",
                    header = T)
NI_A_200.plot <- ggplot(data = NI_A_200.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            size = 0.75) +
  # geom_ribbon(aes(x = years,
  #                 ymin = lcl,
  #                 ymax = ucl,
  #                 group = scenario,
  #                 fill = scenario),
  #             alpha = 0.1) +
  
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

ggsave("figures/NI_A_max_200_1.tiff",
       width = 9,
       height = 6,
       bg = "white") 


# Scenario A all NI 200 pop
NI_A_250.df <- read.csv("data/curlew_HS_NI_A250_max_simple_1.csv",
                        sep = ",",
                        header = T)
NI_A_250.plot <- ggplot(data = NI_A_250.df) +
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

ggsave("figures/NI_A_max_250_ci_1.tiff",
       width = 9,
       height = 6,
       bg = "white") 

# Scenario B all HS_NI_B_means_simple_1
NI_B.df <- read.csv("data/HS_NI_B_means_simple_2.csv",
                    sep = ",",
                    header = T)
NI_B.plot <- ggplot(data = NI_B.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            size = 0.75) +
  # geom_ribbon(aes(x = years,
  #                 ymin = lcl,
  #                 ymax = ucl,
  #                 group = scenario,
  #                 fill = scenario),
  #             alpha = 0.1) +

  
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

ggsave("figures/NI_B_max_150_6.tiff",
       width = 9,
       height = 6,
       bg = "white") 

# Starting pops - 150, 200

# A2 and B2 on same plot, vertical line at 2024 and 2031 because that was where
# looking for stability
# need to extact and do this



#----- RoI  -----
# Scenario A - all - HS_ROI_A_means_simple_1
ROI_A.df <- read.csv("data/HS_ROI_A_means_simple_3.csv",
                    sep = ",",
                    header = T)
ROI.sim.out.plot.ribbon <- ROI.sim.out.plot +
  geom_ribbon(aes(x = years,
                  ymin = lcl,
                  ymax = ucl),
              alpha = 0.3) +

ROI_A.plot <- ggplot(data = ROI_A.df) +
  geom_line(aes(x = years,
                y = pop_size,
                group = scenario,
                colour = scenario),
            size = 0.75) +
  # geom_ribbon(aes(x = years,
  #                 ymin = lcl,
  #                 ymax = ucl,
  #                 group = scenario,
  #                 fill = scenario),
  #             alpha = 0.1) +
  # geom_smooth(aes(x = years,
  #                 y = pop_size,
  #                 group = scenario,
  #                 colour = scenario,
  #                 fill = scenario),
  #             size = 0.25,
  #             span = 0.3) +
  
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

ggsave("figures/ROI_A_max_6.tiff",
       width = 9,
       height = 6,
       bg = "white") 

# Scenario B - all - still need to do this model
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

ggsave("figures/ROI_B_max_ci_2.tiff",
       width = 9,
       height = 6,
       bg = "white") 

# A2 and B2 on same plot, vertical line at 2024 and 2031 because that was where
# looking for stability
# need to extract and do this

NI_A2.df <- read.csv("data/NI_A2_diffpops_clutches_1.csv",
                    sep = ",",
                    header = T)
NI_A2.plot <- ggplot(data = NI_A2.df) +
  facet_wrap(~num_clutches) + 
  geom_line(aes(x = years,
                y = mean/2,
                group = as.factor(initN),
                colour = as.factor(initN)),
            size = 0.75) +
  geom_ribbon(aes(x = years,
                  ymin = lcl/2,
                  ymax = ucl/2,
                  group = as.factor(initN),
                  fill = as.factor(initN)),
              alpha = 0.1) +
  geom_vline(xintercept = 2024,
             linetype = "dashed",
             colour = "grey") +
  geom_vline(xintercept = 2031,
           linetype = "dashed", colour = "grey") +
  
  scale_colour_manual(values = megaPal) +
  scale_fill_manual(values = megaPal) +
  
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

ggsave("figures/NI_B_max_150_5.tiff",
       width = 9,
       height = 6,
       bg = "white") 




#----- Historic Model Plots -----
# NI 
NI.historic.df <- read.csv("data/NI_historic_1987_match750_1.csv",
                           sep = ",",
                           header = TRUE)

pop <- c((5000), (2091), (526))
pop.lcl <- c((3800), (1243), (252))
pop.ucl <- c((6250), (2928), (783))
NI.pop.ests <- data.frame(years = c(1987, 1999, 2013),
                          popsize = pop,
                          pop.min = pop.lcl,
                          pop.max = pop.ucl)

NI.historic.plot <- ggplot(data = NI.historic.df) +
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

ggsave("figures/NI_historic_1.tiff",
       width = 9,
       height = 6,
       bg = "white") 


# RoI 
ROI.historic.df <- read.csv("data/ROI_sims_1985_match400.csv",
                           sep = ",",
                           header = TRUE)

ROIpop <- c((7000), (1000), (138))
ROI.lcl <- c((5000), (743), (100))
ROI.ucl <- c((8500), (1850), (250))
ROI.pop.ests <- data.frame(years = c(1985, 2001, 2015),
                          popsize = ROIpop,
                          pop.min = ROI.lcl,
                          pop.max = ROI.ucl)

ROI.historic.plot <- ggplot(data = ROI.historic.df) +
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

ggsave("figures/ROI_historic_1.tiff",
       width = 9,
       height = 6,
       bg = "white") 



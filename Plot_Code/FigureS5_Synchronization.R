################################################################################
### Figure 6: Oscillation Synchronization ######################################
################################################################################

# oscillations have a period of about 45 ticks

setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path), ".."))

# Oscillate external glutamate to observe synchronization
magnitude = 2
offset = 16

age = 3000

trackers = read.csv(paste0("Final_Pipeline/Final_Out/POST_TRACK/Post_track_",
                           1, ".csv"), header = TRUE)
trackers = trackers[, 2:length(trackers[1, ])]
g_traj = trackers$Out_G_i[c(1990:2389, 2390:2789 + offset)] -
  mean(trackers$Out_G_i[c(1990:2389, 2390:2789 + offset)])

glu_plot_data = data.frame(
  time = 1:800,
  g_traj = g_traj
)

for (repeater_var_10 in 1:4) {
  is.growing = FALSE
  atr = read.csv(paste0("Final_Pipeline/Final_Out/POST_ATR/Post_atr_",
                        repeater_var_10, ".csv"), header = TRUE)
  atr = atr[, 2:length(atr[1, ])]
  atr = as.matrix(atr) # some of my code for ATRs doesn't work on dataframes

  trackers = read.csv(paste0("Final_Pipeline/Final_Out/POST_TRACK/Post_track_",
                             repeater_var_10, ".csv"), header = TRUE)
  trackers = trackers[, 2:length(trackers[1, ])]

  age = length(trackers[, 1])
  unoccupied = which(atr[, 2] == -1)
  occupied = which(is.element(atr[, 1], unoccupied) == FALSE)
  glu_radius = find.glu.radius()

  for (stepwise_timer in (age + 1):(age + 800)) {

    G_m = 30
    G_m_vec = rep(G_m + g_traj[stepwise_timer - age]*magnitude, r*6)

    atr2 = atr

    atr2[, 4:5] = pot.update(atr)[, 4:5]
    atr2[, 6:7] = glu.update(atr)[, 6:7]
    atr2[, 8] = mem.potential.update(atr)[, 8]

    atr = atr2

    atr[unoccupied, c(5, 7, 8)] = 0
    atr[unoccupied, 2:3] = -1

    trackers = update.trackers(stepwise_timer)
    source("BreakScript.R")
  }

  post_traj_glu_osc_1 = trackers$Out_fv_sig[(age + 1):(age + 800)]

  glu_plot_data = cbind(glu_plot_data, post_traj_glu_osc_1[1:800])
  colnames(glu_plot_data)[length(glu_plot_data[1, ])] = paste0("g_osc.",
                                                               repeater_var_10)
}

write.csv(glu_plot_data, paste0("Plot_Code/Plot_Data/Glu_Synch_Data_",
                                magnitude, ".csv"), row.names = FALSE)


# Oscillate external K+ to observe synchronization
# worked with dividing by 3.5 and magnifying by a factor of 1. Now going to
# try other values (change period?) (still effective w smaller magnitude?)
magnitude = 0.05
offset = 25

age = 3000

trackers = read.csv(paste0("Final_Pipeline/Final_Out/POST_TRACK/Post_track_",
                           1, ".csv"), header = TRUE)
trackers = trackers[, 2:length(trackers[1, ])]
k_traj = trackers$Out_K_e[c(2001:2400, 2401:2800 + offset)] -
  mean(trackers$Out_K_e[c(2001:2400, 2401:2800 + offset)])

k_plot_data = data.frame(
  time = 1:800,
  k_traj = k_traj
)

for (repeater_var_10 in 1:4) {
  is.growing = FALSE
  atr = read.csv(paste0("Final_Pipeline/Final_Out/POST_ATR/Post_atr_",
                        repeater_var_10, ".csv"), header = TRUE)
  atr = atr[, 2:length(atr[1, ])]
  atr = as.matrix(atr) # some of my code for ATRs doesn't work on dataframes

  trackers = read.csv(paste0("Final_Pipeline/Final_Out/POST_TRACK/Post_track_",
                             repeater_var_10, ".csv"), header = TRUE)
  trackers = trackers[, 2:length(trackers[1, ])]

  age = length(trackers[, 1])
  unoccupied = which(atr[, 2] == -1)
  occupied = which(is.element(atr[, 1], unoccupied) == FALSE)
  glu_radius = find.glu.radius()

  for (stepwise_timer in (age + 1):(age + 800)) {

    atr[, 4] = atr[, 4] + k_traj[stepwise_timer - age]*magnitude

    atr2 = atr

    atr2[, 4:5] = pot.update(atr)[, 4:5]
    atr2[, 6:7] = glu.update(atr)[, 6:7]
    atr2[, 8] = mem.potential.update(atr)[, 8]

    atr = atr2

    atr[, 4] = atr[, 4] - k_traj[stepwise_timer - age]*magnitude

    atr[unoccupied, c(5, 7, 8)] = 0
    atr[unoccupied, 2:3] = -1

    trackers = update.trackers(stepwise_timer)
    source("BreakScript.R")
  }

  post_traj_k_osc_1 = trackers$Out_fv_sig[(age + 1):(age + 800)]

  k_plot_data = cbind(k_plot_data, post_traj_k_osc_1[1:800])
  colnames(k_plot_data)[length(k_plot_data[1, ])] =
    paste0("k_osc_Gm", G_m, ".", (repeater_var_10 + 1) %% 2 + 1)
}

write.csv(k_plot_data, paste0("Plot_Code/Plot_Data/K_Synch_Data_",
                              magnitude, ".csv"), row.names = FALSE)

k_data = read.csv("Plot_Code/Plot_Data/K_Synch_Data_0.05.csv",
                  header = TRUE)
glu_data = read.csv("Plot_Code/Plot_Data/Glu_Synch_Data_2.csv",
                    header = TRUE)

### Figure S5A: Glutamate synchronization

p1_osc <- glu_data %>% ggplot(aes(x = time)) +
  geom_line(aes(y = rescale(g_traj, c(glu_data[, 3:5])),
                lty = "Scaled Glutamate"), color = "purple", size = 0.8) +
  geom_line(aes(y = g_osc.1, lty = "Signaling Fraction"),
            color = "#F8766D", alpha = 0.7) +
  geom_line(aes(y = g_osc.2, lty = "Signaling Fraction"),
            color = "#F8766D", alpha = 0.7) +
  geom_line(aes(y = g_osc.3, lty = "Signaling Fraction"),
            color = "#F8766D", alpha = 0.7) +
  geom_vline(xintercept = 400.5, color = "slategray", size = 2, alpha = 0.5) +
  scale_linetype_manual(values = c("Signaling Fraction" = "solid",
                                   "Scaled Glutamate" = "21")) +
  labs(title = bquote(bold(A)),
       x = NULL,
       y = "Signaling Fraction",
       linetype = NULL) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11),
    title = element_text(size = 12),
    legend.position = c(0.89, .145),
    legend.margin = margin(0.1, 0, 0, 0.1, unit = "cm"),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.5, "cm")
  )

### Figure S5B: Potassium synchronization

p2_osc <- ggplot(k_data, aes(x = time)) +
  geom_line(aes(y = rescale(k_traj, c(k_data[, 3:5])), lty = "Scaled K+"),
            color = "gold", size = 0.8) +
  geom_line(aes(y = k_osc.1, lty = "Signaling Fraction"),
            color = "#F8766D", alpha = 0.7) +
  geom_line(aes(y = k_osc.2, lty = "Signaling Fraction"),
            color = "#F8766D", alpha = 0.7) +
  geom_line(aes(y = k_osc.3, lty = "Signaling Fraction"),
            color = "#F8766D", alpha = 0.7) +
  geom_vline(xintercept = 400.5, color = "slategray", size = 2, alpha = 0.5) +
  scale_linetype_manual(values = c("Signaling Fraction" = "solid",
                                   "Scaled K+" = "21")) +
  labs(title = bquote(bold(B)),
       x = "Time",
       y = "Signaling Fraction",
       linetype = NULL) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11),
    title = element_text(size = 12),
    legend.position = c(0.89, .145),
    legend.margin = margin(0.1, 0, 0, 0.1, unit = "cm"),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.5, "cm")
  )

p2_osc

p1_osc / p2_osc

combined_plot = p1_osc / p2_osc

ggsave("Plot_Code/Final_Plots/FigureS5_Synchronization.tiff",
       combined_plot, width = 7.5, height = 5, units = "in", dpi = 300)


setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path), ".."))

################################################################################
### Figure5A: Glutamate increase to trigger oscillation depression #############
################################################################################

is.growing = FALSE
atr = read.csv("Output/POST_ATR/Post_atr_14.csv", header = TRUE)
atr = atr[, 2:length(atr[1, ])]
atr = as.matrix(atr) # some of my code for ATRs doesn't work on dataframes

trackers = read.csv("Output/POST_TRACK/Post_track_14.csv", header = TRUE)
trackers = trackers[, 2:length(trackers[1, ])]

r = 150 # radius of network
G_m = 30 # basal glutamate concentration in media
G_m_vec = rep(G_m, r*6 + 6)

age = length(trackers[, 1])
unoccupied = which(atr[, 2] == -1)
occupied = which(is.element(atr[, 1], unoccupied) == FALSE)
glu_radius = find.glu.radius()

for (stepwise_timer in (age + 1):(age + 400)) {
  if (stepwise_timer < (age + 300) & stepwise_timer > (age + 100)) {
    G_m = 35
    G_m_vec = rep(G_m, r*6 + 6)
  } else {
    G_m = 30
    G_m_vec = rep(G_m, r*6 + 6)
  }

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

post_traj_glu_increase_1 = trackers$Out_fv_sig[(age + 1):(age + 400)]

write.csv(trackers, "Plot_Code/Plot_Data/Glu_Increase_Trackers.csv",
          row.names = FALSE)

trackers = read.csv("Plot_Code/Plot_Data/Glu_Increase_Trackers.csv",
                    header = TRUE)

p1 = trackers[(age + 1):(age + 400), ] %>% ggplot(aes(x = Time - age)) +
  geom_vline(xintercept = 100, color = "purple", size = 1, lty = "dashed") +
  geom_vline(xintercept = 300, color = "purple", size = 1, lty = "dashed") +
  geom_line(aes(y = Tot_fv_sig), color = "#F8766D", size = 0.4) +
  labs(title = bquote(bold(A) ~ "  Increased glutamate"),
       x = "",
       y = "Total Signaling Fraction") +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)
  )

################################################################################
### Figure5B: K+ increased to trigger oscillations #############################
################################################################################

trackers2 = read.csv("Output/POST_TRACK/Post_track_KShock_750.csv",
                     header = TRUE)

k_shock_begin = 750

trackers_trimmed = trackers2[(k_shock_begin - 200):(k_shock_begin + 800), ]

p2 = ggplot() +
  geom_line(aes(x = trackers_trimmed$Time - k_shock_begin + 200,
                y = trackers_trimmed$Tot_fv_sig),
            color = "#F8766D", size = 0.4) +
  geom_rect(aes(xmin = 200, xmax = 204, ymin = 0, ymax = Inf),
            #color = "slategray", fill = "slategray", alpha = 0.5) +
            color = "gold", fill = "gold", alpha = 0.5) +
  labs(title = bquote(bold(B) ~ "  Increased K+"),
       x = "",
       y = "") +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)
  )

################################################################################
### FIGURE 5C GROWTH WITH LOW GLUTAMATE ########################################
################################################################################

# lg = low glutamate
trackers_lg = read.csv("Output/POST_TRACK/Post_track_LowGlu.csv", header = TRUE)

trackers_lg = trackers_lg %>% filter(Time > 780 & Time <= 1580)

# hg = high glutamate
trackers_hg = read.csv("Output/POST_TRACK/Post_track_1.csv", header = TRUE)
trackers_hg = trackers_hg %>% filter(Time <= 1510 & Time > 710)

scale_factor = max(trackers_hg$Radius)/max(trackers_lg$Out_fv_sig)

trackers_combined = data.frame(
  Time = 1:800,
  Out_fv_sig_lg = trackers_lg$Tot_fv_sig,
  Out_fv_sig_hg = trackers_hg$Tot_fv_sig,
  Radius_lg = trackers_lg$Radius,
  Radius_hg = trackers_hg$Radius
)

linetypes = c("Normal" = "solid", "Reduced" = "11")

p3 = ggplot(trackers_combined, aes(x = Time)) +
  geom_line(aes(y = Out_fv_sig_hg*max(Radius_hg)/max(Out_fv_sig_lg),
                lty = "Normal"),
            color = "#F8766D", size = 0.4) +
  geom_line(aes(y = Radius_hg, lty = "Normal"), color = "slategray",
            size = 0.5) +
  geom_line(aes(y = Out_fv_sig_lg*max(Radius_hg)/max(Out_fv_sig_lg),
                lty = "Reduced"),
            color = "#F8766D", size = 0.4) +
  geom_line(aes(y = Radius_lg, lty = "Reduced"), color = "slategray",
            size = 0.5) +
  scale_linetype_manual(values = linetypes) +
  scale_y_continuous(
    name = "Radius (Cells)",
    sec.axis = sec_axis(trans = ~./scale_factor,
                        name = "Total Signaling Fraction")
  ) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y.left = element_text(color = "slategray"),
    axis.text.y.left = element_text(color = "slategray"),
    axis.title.y.right = element_text(color = "#F8766D"),
    axis.text.y.right = element_text(color = "#F8766D"),
    text = element_text(size = 12),
    legend.position = c(0.92, 0.15)
  ) +
  labs(
    x = "Time",
    title = bquote(bold(C) ~ "  Reduced glutamate"),
    linetype = NULL
  )

### Combined plot

p4 = grid.arrange(
  arrangeGrob(
    p1, p2, ncol = 2
  ),
  p3, ncol = 1
)

ggsave("Plot_Code/Final_Plots/Figure5_MediaPerturbations.tiff", p4,
       width = 7.5, height = 6.5, units = "in", dpi = 300)


################################################################################
### FIGURE: EFFECTIVENESS OF SIGNALING #########################################
################################################################################

warning("Must be run after 'Figure7_SignalingEfficacy.R' to prevent errors")

setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path), ".."))

mg = read.csv(
  "Final_Pipeline/Final_Out/Other/Minimum_Glu_Signaling_Thresh[0,1]_SD3.csv",
  header = TRUE)
colnames(mg) = mg[1, ]
mg = mg[2:length(mg[, 1]), ]
mean_glu = colMeans(mg)
vals = rep(0, length(atr[, 1]))
vals[as.numeric(colnames(mg))] = mean_glu
write.csv(vals, "Plot_Code/Plot_Data/Efficacy_SigVals_013.csv",
          row.names = FALSE)

mg = read.csv(
  "Final_Pipeline/Final_Out/Other/Minimum_Glu_Signaling_Thresh[1,3]_SD3.csv",
  header = TRUE)
colnames(mg) = mg[1, ]
mg = mg[2:length(mg[, 1]), ]
mean_glu = colMeans(mg)
vals = rep(0, length(atr[, 1]))
vals[as.numeric(colnames(mg))] = mean_glu
write.csv(vals, "Plot_Code/Plot_Data/Efficacy_SigVals_133.csv",
          row.names = FALSE)


################################################################################
################################################################################
################################################################################

vals_000 = data.frame(
  "ID" = 1:length(atr[, 1]),
  "Value" = c(read.csv("Plot_Code/Plot_Data/Efficacy_NonVals.csv",
                       header = TRUE))[[1]])
vals_031 = data.frame(
  "ID" = 1:length(atr[, 1]),
  "Value" = c(read.csv("Plot_Code/Plot_Data/Efficacy_SigVals.csv",
                       header = TRUE))[[1]])
vals_013 = data.frame(
  "ID" = 1:length(atr[, 1]),
  "Value" = c(read.csv("Plot_Code/Plot_Data/Efficacy_SigVals_013.csv",
                       header = TRUE))[[1]])
vals_133 = data.frame(
  "ID" = 1:length(atr[, 1]),
  "Value" = c(read.csv("Plot_Code/Plot_Data/Efficacy_SigVals_133.csv",
                       header = TRUE))[[1]])


traj_000 = read.csv(
  "Final_Pipeline/Final_Out/POST_TRACK/Post_track_Thresh[0,0]_SD0.csv",
  header = TRUE)
traj_031 = read.csv(
  "Final_Pipeline/Final_Out/POST_TRACK/Post_track_1.csv",
  header = TRUE)
traj_013 = read.csv(
  "Final_Pipeline/Final_Out/POST_TRACK/Post_track_Thresh[0,1]_SD3.csv",
  header = TRUE)
traj_133 = read.csv(
  "Final_Pipeline/Final_Out/POST_TRACK/Post_track_Thresh[1,3]_SD3.csv",
  header = TRUE)

p1a = ggplot(vals_000 %>% filter(Value > 0)) +
  geom_density(aes(x = Value), color = "slategray", fill = "purple",
               alpha = 0.5) +
  annotate("text", x = 0,
           y = max(density((vals_000 %>% filter(Value > 0))$Value)$y)*0.9,
           label = "Threshold = 0",
           color = "black", size = 3.5, hjust = 0) +
  labs(
    x = "Time-Averaged Internal Glutamate",
    y = "Cell Frequency",
    title = bquote(bold(A)) # Threshold = 0"
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11)
  )

p2a = ggplot(vals_013 %>% filter(Value > 0)) +
  geom_density(aes(x = Value), color = "slategray", fill = "purple",
               alpha = 0.5) +
  annotate("text", x = 0,
           y = max(density((vals_013 %>% filter(Value > 0))$Value)$y)*0.9,
           label = "Threshold ~ [0, 1]",
           color = "black", size = 3.5, hjust = 0) +
  labs(
    x = "Time-Averaged Internal Glutamate",
    y = "Cell Frequency",
    title = bquote(bold(C)) # Threshold = [0, 1]"
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11)
  )

p3a = ggplot(vals_031 %>% filter(Value > 0)) +
  geom_density(aes(x = Value), color = "slategray", fill = "purple",
               alpha = 0.5) +
  annotate("text", x = 0,
           y = max(density((vals_031 %>% filter(Value > 0))$Value)$y)*0.9,
           label = "Threshold ~ [0, 3]",
           color = "black", size = 3.5, hjust = 0) +
  labs(
    x = "Time-Averaged Internal Glutamate",
    y = "Cell Frequency",
    title = bquote(bold(E)) # Threshold = [0, 3]"
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11)
  )

p4a = ggplot(vals_133 %>% filter(Value > 0)) +
  geom_density(aes(x = Value), color = "slategray", fill = "purple",
               alpha = 0.5) +
  annotate("text", x = 0,
           y = max(density((vals_133 %>% filter(Value > 0))$Value)$y)*0.9,
           label = "Threshold ~ [1, 3]",
           color = "black", size = 3.5, hjust = 0) +
  labs(
    x = "Time-Averaged Internal Glutamate",
    y = "Cell Frequency",
    title = bquote(bold(G)) # Threshold = [1, 3]"
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11)
  )

################################################################################

p1b = ggplot(traj_000[2001:2200, ]) +
  geom_line(aes(x = Time - 2000, y = Tot_fv_sig), color = "#F8766D") +
  labs(
    x = "Time",
    y = "Signaling Fraction",
    title = bquote(bold(B))
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11)
  )

p2b = ggplot(traj_013[2001:2200, ]) +
  geom_line(aes(x = Time - 2000, y = Tot_fv_sig), color = "#F8766D") +
  labs(
    x = "Time",
    y = "Signaling Fraction",
    title = bquote(bold(D))
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11)
  )

p3b = ggplot(traj_031[2001:2200, ]) +
  geom_line(aes(x = Time - 2000, y = Tot_fv_sig), color = "#F8766D") +
  labs(
    x = "Time",
    y = "Signaling Fraction",
    title = bquote(bold(F))
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11)
  )

p4b = ggplot(traj_133[2001:2200, ]) +
  geom_line(aes(x = Time - 2000, y = Tot_fv_sig), color = "#F8766D") +
  labs(
    x = "Time",
    y = "Signaling Fraction",
    title = bquote(bold(H))
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11)
  )

# Create titles for each row
row_titles <- list(
  ggplot() + ggtitle("Threshold = 0") + theme_void(),
  ggplot() + ggtitle("Threshold = [0, 1]") + theme_void(),
  ggplot() + ggtitle("Threshold = [0, 3]") + theme_void(),
  ggplot() + ggtitle("Threshold = [1, 3]") + theme_void()
)

final_plot =
  (p1a + plot_spacer() + p1b) /
  (p2a + plot_spacer() + p2b) /
  (p3a + plot_spacer() + p3b) /
  (p4a + plot_spacer() + p4b)

final_plot = p1a + plot_spacer() + p1b + p2a + plot_spacer() + p2b +
  p3a + plot_spacer() + p3b + p4a + plot_spacer() + p4b +
  plot_layout(nrow = 4, ncol = 3, widths = c(1, 0.03, 1))

ggsave("Plot_Code/Final_Plots/FigureS7_VarThresholds.tiff", final_plot,
       width = 7.5, height = 7.5, units = "in", dpi = 300)



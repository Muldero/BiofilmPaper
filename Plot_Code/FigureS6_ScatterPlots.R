setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path), ".."))

sig_traj = readRDS(
  "Output/Other/Membrane_Potential_Matrix_1.RDS"
)[1501:3000, ]

atr = read.csv("Output/POST_ATR/Post_atr_1.csv",
               header = TRUE)[, 2:18]

external = readRDS(
  "Output/Other/Membrane_Potential_External_1.RDS")
occupied = which(colMeans(sig_traj) != 0)

ex_occ = intersect(external, occupied)

glu_vals = read.csv("Plot_Code/Plot_Data/Efficacy_SigVals.csv", header = TRUE)
glu_vals = c(glu_vals)[[1]]


# identify signaling cells
for (i in 1:1500) {
  atr_temp = atr
  atr_temp[, 8] = sig_traj[i, ]

  sig_traj[i, ] = 0
  sig_traj[i, find.sigs.fv(at_tab = atr_temp)] = 1
}

data_ex = data.frame(
  Ids = ex_occ,
  Sig_Rate = colMeans(sig_traj[, ex_occ]),
  Threshold = atr[ex_occ, 2],
  Glu = glu_vals[ex_occ],
  Depth = atr[ex_occ, 9]
)

data_occ = data.frame(
  Ids = occupied,
  Sig_Rate = colMeans(sig_traj[, occupied]),
  Threshold = atr[occupied, 2],
  Glu = glu_vals[occupied],
  Depth = max(atr[occupied, 9]) - atr[occupied, 9]
)

data_occ$Depth[which.max(data_occ$Depth)[1]] = 140

### Figure S6A: Internal glutamate x Depth

p1 = ggplot(data_occ) +
  geom_point(mapping = aes(x = Depth, y = Glu, color = Sig_Rate),
             size = 0.5, alpha = 0.8) +
  scale_color_gradientn(
    colors = c("black", "#00CCCC", "cyan"), # 90% cyan at 0.5, full cyan at 1
    values = c(0, 0.5, 1),                 # Breakpoints
    limits = c(0, 1),                       # Ensure scale spans 0 to 1
    breaks = c(0, 0.5, 1)              # Set legend breakpoints at 0 and 1
  ) +
  scale_x_reverse() +
  labs(title = bquote(bold(A)),
       x = "Depth within the Biofilm",
       y = "Internal Glutamate",
       color = "% Time \n Signaling") +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    legend.key.size = unit(0.4, "cm"),         # Legend key size
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = c(0.99, 0.01),
    legend.justification = c(1, 0)
  )

### Figure S6B: Internal glutamate x Stress threshold

p2 = ggplot(data_occ) +
  geom_point(mapping = aes(x = Threshold, y = Glu, color = Depth),
             size = 0.5) +
  scale_color_gradient(
    low = "lightgray",
    high = "black",
    trans = "reverse",
    breaks = c(0, 70, 140)
  ) +
  labs(title = bquote(bold(B)),
       x = "Stress Threshold",
       y = "Internal Glutamate") +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    legend.key.size = unit(0.4, "cm"),         # Legend key size
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = c(0.99, 0.01),
    legend.justification = c(1, 0)
  )

### Figure S6C: Internal glutamate x Time signaling

p3 = ggplot(data_occ) +
  geom_point(mapping = aes(x = Sig_Rate, y = Glu, color = Depth),
             size = 0.5) +
  scale_color_gradient(low = "lightgray", high = "black", trans = "reverse") +
  labs(title = bquote(bold(C)),
       x = "% Time Signaling",
       y = "Internal Glutamate") +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    legend.key.size = unit(0.4, "cm"),         # Legend key size
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = "none",
    legend.justification = c(1, 0)
  )

### Figure S6D: Time signaling x Stress threshold

p4 = ggplot(data_occ) +
  geom_point(mapping = aes(x = Threshold, y = Sig_Rate, color = Depth),
             size = 0.5) +
  scale_color_gradient(low = "lightgray", high = "black", trans = "reverse") +
  labs(title = bquote(bold(D)),
       x = "Stress Threshold",
       y = "% Time Signaling") +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    legend.key.size = unit(0.4, "cm"),         # Legend key size
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.position = "none",
    legend.justification = c(1, 0)
  )

combined = (p1 + p2)/(p3 + p4)

ggsave(filename = "Plot_Code/Final_Plots/FigureS6_ScatterPlots.tiff",
       plot = combined, width = 7.5, height = 6.5, units = "in", dpi = 300)

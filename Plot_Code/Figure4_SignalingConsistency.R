etwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path), ".."))

frac_osc_sig = c()

for (i in 1:5) {
  sig_traj = readRDS(paste0(
    "Output/Other/Membrane_Potential_Matrix_", i, ".RDS")
    )[1501:3000, ]

  trackers = readRDS(paste0(
    "Output/Other/Membrane_Potential_Tracker_", i, ".RDS"))

  external = readRDS(paste0(
    "Output/Other/Membrane_Potential_External_", i, ".RDS"))
  occupied = which(colMeans(sig_traj) != 0)

  trackers = filter(trackers, trackers$Time > (3000 - length(sig_traj[, 1])))
  peaks = get.only.max.peaks(trackers$Out_fv_sig)
  dips = get.only.max.peaks(trackers$Out_fv_sig, midline = 0.25, top = FALSE)
  num_osc = length(dips)

  temp_store = 1:length(external)

  for (j in 1:length(external)) {
    k = external[j]
    traj_peaks = which(sig_traj[, k] < -390)

    osc_count = 0

    for (x in 1:(num_osc - 1)) {
      if (sum(traj_peaks > dips[x] & traj_peaks < dips[x + 1]) > 0) {
        osc_count = osc_count + 1
      }
    }

    temp_store[j] = osc_count/(num_osc - 1)

    if (temp_store[j] > 1) {
      temp_store[j] = 1
    }
  }

  frac_osc_sig = c(frac_osc_sig, temp_store)
}

fos_df = data.frame(frac = frac_osc_sig)
hp = ggplot(fos_df, aes(x = frac)) +
  geom_histogram(breaks = seq(0, 1, length.out = 12),
                 fill = "#F8766D", color = "red", alpha = 0.9) +
  labs(title = "",
       x = "Fraction of Oscillations Spent Signaling",
       y = "Density") +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    axis.text.y = element_blank()
  )


ggsave(filename = "Plot_Code/Final_Plots/Figure4_SignalingConsistency.tiff",
       plot = hp, width = 5.2, height = 4, units = "in", dpi = 300)

################################################################################
### FIGURE 2A and C ############################################################
################################################################################

# Figure 2 B and D were generated elsewhere and the figure was assembled in
# powerpoint

setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path), ".."))

### Figure 2A: In vitro traces

traces = read.csv("Plot_Code/Plot_Data/InVitroTraces.csv", header = TRUE)[, 2:5]
colnames(traces) = c("Time", "Inner", "Outer", "Total")
traces = traces %>% filter(Time > 7)

normalizer = function(vector) {
  vector = vector - min(vector)
  return(vector/max(vector))
}

traces = traces %>% filter(Time > 9.5, Time < 27.6)

range_out_iv = max(traces$Outer) - min(traces$Outer)
range_in_iv = max(traces$Inner) - min(traces$Inner)

p1_IV = ggplot(traces %>% filter(Time > 9.5, Time < 27.6),
               aes(x = Time - 9.5)) +
  geom_line(aes(y = normalizer(normalizer(Inner) - Time/28)*
                  range_out_iv + min(Outer), lty = "Inner"),
            color = "cyan", size = 1) +
  geom_line(aes(y = Outer, lty = "Outer"), color = "cyan", size = 1) +
  scale_y_continuous(
    name = "Outer ThT Intensity",
    sec.axis = sec_axis(trans = ~((. - min(traces$Outer))/range_out_iv)*range_in_iv +
                          min(traces$Inner),
                        name = "Adjusted Inner ThT Intensity")
  ) +
  labs(
    x = "Time (hrs)",
    y = "Scaled ThT Intensity",
    title = bquote(bold("A")),
    linetype = bquote(italic("In vitro"))
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 10),
    legend.position = c(.905, .27),
    legend.spacing = unit(0, "cm"),
    legend.margin = margin(0.1, 0, 0, 0.1, unit = "cm"),
    legend.text = element_text(size = 8)
  )

### Figure 2C: Simulated traces

trackers = read.csv("Final_Pipeline/Final_Out/POST_TRACK/Post_track_2.csv",
                    header = TRUE)

trackers = trackers %>% filter(Time < 2730, Time > 2500)

range_out = max(trackers$Out_fv_sig) - min(trackers$Out_fv_sig)
range_in = max(trackers$In_fv_sig) - min(trackers$In_fv_sig)

p2_SIM = ggplot(trackers, aes(x = Time - 2500)) +
  geom_line(aes(y = normalizer(In_fv_sig)*range_out + min(Out_fv_sig),
                lty = "Inner"), color = "#F8766D",
            size = 1) +
  geom_line(aes(y = Out_fv_sig, lty = "Outer"), color = "#F8766D",
            size = 1) +
  scale_y_continuous(
    name = "Outer Fraction Signaling",
    sec.axis = sec_axis(trans = ~((. - min(trackers$Out_fv_sig))/
                                    range_out)*range_in +
                          min(trackers$In_fv_sig),
                        name = "Inner Fraction Signaling")
  ) +
  labs(
    x = "Time (ticks)",
    y = "Scaled Fraction Signaling",
    title = bquote(bold("C")),
    linetype = bquote("Simulated")
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 10),
    legend.position = c(.905, .27),
    legend.spacing = unit(0, "cm"),
    legend.margin = margin(0.1, 0, 0, 0.1, unit = "cm"),
    legend.text = element_text(size = 8)
  )

combined = p1_IV/p2_SIM

combined

ggsave("Plot_Code/Final_Plots/Figure2AC_InVistroSimulated.tiff", combined,
       width = 5.5, height = 4.5, units = "in", dpi = 300)


################################################################################
### FIGURE 3: GROWTH ILLUSTRATION ##############################################
################################################################################

setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path), ".."))

trackers = read.csv("Output/POST_TRACK/Post_track_1.csv",
                    header = TRUE) %>% filter(Time < 1599 & Time > 600)

scale_factor = max(trackers$Radius)/max(trackers$Out_fv_sig)

p1_grow = ggplot(trackers, aes(x = Time - 600)) +
  geom_line(aes(y = Out_fv_sig*scale_factor), color = "#F8766D",
            lty = "solid", size = 0.3) +
  geom_line(aes(y = Radius), color = "slategray",
            lty = "solid", size = 0.3) +
  scale_y_continuous(
    name = "Radius (Cells)",
    sec.axis = sec_axis(trans = ~./scale_factor,
                        name = "Outer Signaling Fraction")
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y.left = element_text(color = "slategray"),
    axis.text.y.left = element_text(color = "slategray"),
    axis.title.y.right = element_text(color = "#F8766D"),
    axis.text.y.right = element_text(color = "#F8766D"),
    text = element_text(size = 12)
  ) +
  labs(x = "Time")

ggsave("Plot_Code/Final_Plots/Figure3_Growth.tiff", p1_grow,
       width = 5.2, height = 3, units = "in", dpi = 300)


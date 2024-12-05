atr = read.csv("Output/POST_ATR/Post_atr_16.csv", header = TRUE)
atr = atr[, 2:length(atr[1, ])]
atr = as.matrix(atr)

atr1 = atr

F_v = Eff(atr1[which(atr1[, 2] != -1), 8])
mid_val = weighted.mean(c(mean(F_v), 1), c(2/3, 1/3))

# Now identify the maxima for cells below the midpoint
lower_half_fv = F_v[which(F_v < mid_val)]
dens = density(lower_half_fv, adjust = 0.2)
xval = dens$x[which.max(dens$y)]

# Now define our cutoff to be a weighted average between x and midval
cutoff = weighted.mean(c(xval, 1), c(3/4, 1/4))

F_v = F_v %>% as.data.frame()
colnames(F_v) = "V"

dat = with(density(F_v$V, from = min(F_v), to = max(F_v)), data.frame(x, y))

p_mem = ggplot(data = dat, mapping = aes(x = x, y = y)) +
  geom_vline(aes(xintercept = cutoff), size = 1, linetype = "dotted") +
  geom_area(mapping = aes(x = ifelse(x <= cutoff, x, 0),
                          fill = "Nonsignaling"), alpha = 0.9) +
  geom_area(mapping = aes(x = ifelse(x >= cutoff, x, 0), fill = "Signaling"),
            alpha = 0.9) +
  geom_line(color = "black") +
  scale_fill_manual(values = c("Signaling" = "slategray",
                                "Nonsignaling" = "lightgray")) +
  xlim(min(F_v), max(F_v)) +
  labs(title = "",
       x = "F(V)",
       y = "Density",
       fill = NULL) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    legend.position = c(0.86, 0.877)
  )

ggsave("Plot_Code/Final_Plots/FigureS2_SignalingDefinition.tiff", p_mem,
       width = 5.2, height = 4, units = "in",
       dpi = 300)

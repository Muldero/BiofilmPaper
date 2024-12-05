################################################################################
### FIGURE7: EFFECTIVENESS OF SIGNALING ########################################
################################################################################

# Note the B and D are not printed in this plot and were added in post

setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path), ".."))

mg_nosig = read.csv(
  "Final_Pipeline/Final_Out/Other/Minimum_Glu_Signaling_Thresh[0,0]_SD0.csv",
  header = TRUE)
mg_sig = read.csv(
  "Final_Pipeline/Final_Out/Other/Minimum_Glu_Signaling_1.csv",
  header = TRUE)

colnames(mg_nosig) = mg_nosig[1, ]
colnames(mg_sig) = mg_sig[1, ]


mg_nosig = mg_nosig[2:length(mg_nosig[, 1]), ]
mg_sig = mg_sig[2:length(mg_sig[, 1]), ]

mean_glu_s = colMeans(mg_sig)
mean_glu_ns = colMeans(mg_nosig)

nonsig_vals = rep(0, length(atr[, 1]))
nonsig_vals[as.numeric(colnames(mg_nosig))] = mean_glu_ns#*10

sig_vals = rep(0, length(atr[, 1]))
sig_vals[as.numeric(colnames(mg_sig))] = mean_glu_s#*10

write.csv(sig_vals, "Plot_Code/Plot_Data/Efficacy_SigVals.csv",
          row.names = FALSE)
write.csv(nonsig_vals, "Plot_Code/Plot_Data/Efficacy_NonVals.csv",
          row.names = FALSE)


nonsig_vals = read.csv("Plot_Code/Plot_Data/Efficacy_NonVals.csv",
                       header = TRUE)
nonsig_vals = c(nonsig_vals)[[1]]
sig_vals = read.csv("Plot_Code/Plot_Data/Efficacy_SigVals.csv", header = TRUE)
sig_vals = c(sig_vals)[[1]]

nonsig_df = as.data.frame(cbind(1:length(atr[, 1]), nonsig_vals))
colnames(nonsig_df) = c("Id", "Value")

sig_df = as.data.frame(cbind(1:length(atr[, 1]), sig_vals))
colnames(sig_df) = c("Id", "Value")

purple_palette = colorRampPalette(c("#FFFFFF", "purple"))(100)

# Figure 7A: Density of internal glutamate in nonsignaling biofilm

p1 = ggplot(nonsig_df %>% filter(Value > 0)) +
  geom_density(aes(x = Value), color = "slategray", fill = "purple",
               alpha = 0.5) +
  labs(
    x = "Time-Averaged Internal Glutamate",
    y = "Cell Frequency",
    title = bquote(bold(A) ~ "  Without Signaling")
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)
  )

# Figure 7B: Glutamate distribution across nonsignaling biofilm

p2 = plot.hex.net.gg(cat = nonsig_vals*10, altcols = purple_palette,
                     are.unoccupied = which(nonsig_vals == 0), remake = TRUE)

# Figure 7C: Density of internal glutamate in signaling biofilm

p3 = ggplot(sig_df %>% filter(Value > 0)) +
  geom_density(aes(x = Value), color = "slategray", fill = "purple",
               alpha = 0.5) +
  labs(
    x = "Time-Averaged Internal Glutamate",
    y = "Cell Frequency",
    title = bquote(bold(C) ~ "  With Signaling")
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)
  )

# Figure 7D: Glutamate distribution across signaling biofilm

p4 = plot.hex.net.gg(cat = sig_vals*10, altcols = purple_palette,
                     are.unoccupied = which(sig_vals == 0), remake = TRUE)

# Empty plot for padding
empty_plot = ggplot() + theme_void()

# Arrange plots with empty plot for padding
p5 = grid.arrange(p1, p2, empty_plot, empty_plot, p3, p4,
                  nrow = 3,
                  heights = c(1, 0.1, 1))

ggsave("Plot_Code/Final_Plots/Figure7_SignalingEfficacy.tiff", p5,
       width = 7.5, height = 7.5, units = "in", dpi = 600)



### Calculate Gini coefficient

library(ineq)

Gini(nonsig_vals)
Gini(sig_vals)


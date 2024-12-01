# This script is sourced in loops where the "stop" button in RStudio fails.
# By uncommenting "break()" here the loop will stop.



# if (length(find.sigs.fv())/length(occupied) > 0.42 & is.growing == FALSE) {
#   saveRDS(atr, "Final_Pipeline/Final_Out/Peaks_ATR7.RDS")
#   plot.hex.net.quickndirty()
#   print("done")
  # break()
# }

# print(summary(atr[occupied, 2]))

# is.growing = TRUE
do.plot = FALSE

# print(max(atr[occupied, 9]))

# do.plot = TRUE
#
# max.smoothing = function(v, d) {
#   v2 = v
#   for (i in (d + 1):(length(v) - d)) {
#     r_temp = v[(i - d):(i + d)]
#     m1 = max(r_temp)
#     m2 = max(r_temp[-which(r_temp == m1)[1]])
#     v2[i] = mean(m1, m2)
#   }
#   return(v2)
# }

# plot.singlecell.trajs()

# plot.hex.net.quickndirty()

# source("Functions_Hunger3.8_SimulUpdate.R", verbose = FALSE)

# break()

# print(stepwise_timer)
# print(summary(atr[, 7]))

# plot(density(atr[, 7]), xlim = c(0, 10))

# print(summary(atr[occupied, 7]))

# par(mfrow = c(1, 1))
# if (stepwise_timer %% 20 == 0) {
#   plot(trackers$Out_fv_sig[3001:stepwise_timer], type = "l", ylim = c(0, 0.65))
#   # line(pot_traj[1:(stepwise_timer - 3000)] + 0.25)
#   abline(h = 0.43)
# }


# if (stepwise_timer %in% c(3300:3400)) {
#   plot.a.tiff.please(stepwise_timer, plot.boundary = TRUE)
# }



# plot(trackers$Out_sig[1:stepwise_timer], type = "l")
# # plot(trackers$Out_sig_inher[1:stepwise_timer], type = "l")
# # plot(trackers$Out_non_inher[1:stepwise_timer], type = "l")
# plot(density(Eff(atr[find.external(), 8])))
# F_v = Eff(atr[occupied, 8])
# mid_val = mean(c(mean(F_v), 1))
# abline(v = mid_val)
#
# par(mfrow = c(1, 1))

# print(trackers$Out_fv_sig[stepwise_timer])

# print(plot.hex.net.gg())

# plot(density(atr[occupied, 7]))

# par(mfrow = c(2, 1))
# plot(trackers$Out_fv_sig[1:stepwise_timer], type = "l")
# F_v = Eff(atr[occupied, 8])
# plot(density(Eff(atr[occupied, 8])))
#
# F_v = Eff(atr[occupied, 8])
# dens = density(F_v)
# local_maxima = which(diff(sign(diff(dens$y))) < 0) + 1
# mode = sort(dens$y[local_maxima], decreasing = TRUE)[1]
# midpoint = dens$x[which(is.element(dens$y, mode))]
# mid_val = (midpoint + max(F_v))/2
#
# abline(v = mid_val)
# # abline(v = mean(F_v))
# # abline(v = Eff(-160), lty = "dotted")
# par(mfrow = c(1, 1))

# plot.singlecell.trajs()

# print(mean(atr[occupied, 4]))

# plot.hex.net.quickndirty(double = FALSE, remake = FALSE)

# plot(atr[, 9], atr[, 4])

# print(summary(atr[occupied, 8]))
# print(summary(atr[occupied, 7]))
#
# print(sum(atr[, 7] < atr[, 2]))
# print(mean(Eff(atr[occupied, 8])))
#

# source("Functions_Hunger3.6.R")

# if (osc_start == TRUE) {
#   prop_occupied = 0.001
# break()
# }

# print(atr[1, 7])

# if (stepwise_timer %% 10 == 0) {
# print(paste0("Tick: ", stepwise_timer, "  |  Size: ", max(atr[occupied, 9])))
# #
# #   Sys.sleep(2)
# #   flexy.plot(atr[, 7])
# #   Sys.sleep(10)
# }

# print(trackers[stepwise_timer, ]$Out_fv_sig)

# print(warnings())

print(paste0(stepwise_timer))

# print(stepwise_timer)

### Let's do this! #############################################################
# atr = attribute table, records all information about the biofilm
# at a given tick

# trackers = records variety of data from atr for each tick
# across the entire run
################################################################################

# Initialize the biofilm

parent_matrix = cbind(1:7, c(NA, 1, 1, 1, 1, 1, 1)) # parents for cells in atr
colnames(parent_matrix) = c("Offspring", "Parent")

if (grow.all.the.way == TRUE) {
  atr = grow.a.biofilm2(method = "entire")
} else {
  atr = grow.a.biofilm2(method = "first_layer") # make.hex.atr()
  for (i in 1:ceiling(init_radius/gens_per_tick)) {
    atr = grow.a.biofilm2(method = "add_layer")
  }
}
################################################################################

# rows of the atr representing empty hexes
unoccupied = which(atr[, 2] == -1)

# rows of the atr with cells in them
occupied = which(is.element(atr[, 1], unoccupied) == FALSE)

# boundary between interior and exterior
glu_radius = find.glu.radius()

# record potassium, glutamate and memebrane potential for 100 random cells
pot_traj = matrix(nrow = 100, ncol = 0)
glu_traj = matrix(nrow = 100, ncol = 0)
mem_pot_traj = matrix(nrow = 100, ncol = 0)

traj.created = FALSE
is.growing = TRUE

# record signaling state of every cell in the biofilm (memory intensive!)
if (track.all.sigs == TRUE) {
  sig_tracking = matrix(0, nrow = ticks, ncol = length(atr[, 1]))
}

# Now for the actual simulation
for (stepwise_timer in 1:ticks) {

  # Do we do a potassium shock?
  if (stepwise_timer >= min(K_shock) & stepwise_timer < max(K_shock + 5)) {
    K_m = kshock_level
  } else {
    K_m = 8
  }

  # Simultaneously update potassium, glutamate, and membrane potential
  atr2 = atr

  atr2[, 4:5] = pot.update(atr)[, 4:5]
  atr2[, 6:7] = glu.update(atr)[, 6:7]
  atr2[, 8] = mem.potential.update(atr)[, 8]

  atr = atr2

  # Correct values in unoccupied cells
  atr[unoccupied, c(5, 7, 8)] = 0
  atr[unoccupied, 2:3] = -1


  # While growing
  if (length(occupied)/length(atr[, 1]) < prop_occupied) {
    # growth the biofilm by 1/40 of a layer
    atr = grow.a.biofilm2(method = "add_layer")

    # update occupied cells and inner/outer boundary
    unoccupied = which(atr[, 2] == -1)
    occupied = which(is.element(atr[, 1], unoccupied) == FALSE)
    glu_radius = find.glu.radius()

  } else if (is.growing == TRUE) {
    # first tick after growth stops perform some stuff

    is.growing = FALSE # stop growth

    # make signal recording matrix
    who_is_signaling = matrix(nrow = length(occupied), ncol = 1)
    colnames(who_is_signaling) = 1

    # make signaling trajectory matrix
    traj.created = TRUE
    traj_time_init = stepwise_timer

    if (track.all.sigs == FALSE) {
      signal_trajectories = matrix(0, nrow = ticks - stepwise_timer + 1,
                                   ncol = 1000)
    }

    cells_to_track = sample(find.external(), 1000)

    if (track.min.glu == TRUE) {
      # make minimum glutamate tracking matrix
      min_glu = atr[occupied, 7]
    }

    # all of the above are just various recording mechanisms
  } else {
    # after growth starts

    # update signaling tracker
    who_is_signaling = check.whos.signaling()

    if (track.all.sigs == FALSE) {
      # update trajectory tracker
      signal_trajectories[stepwise_timer - traj_time_init + 1, ] =
        atr[cells_to_track, 8]
    }

    if (track.min.glu == TRUE) {
      # update minimum glutamate tracker
      min_glu = rbind(min_glu, atr[occupied, 7])
    }

    # Update single cell traces
    pot_traj = cbind(pot_traj, atr[cells_to_track[1:100], 5])
    glu_traj = cbind(glu_traj, atr[cells_to_track[1:100], 7])
    mem_pot_traj = cbind(mem_pot_traj, atr[cells_to_track[1:100], 8])
  }

  if (track.all.sigs == TRUE) {
    sig_tracking[stepwise_timer, ] = atr[, 8]
  }

  # update trackers (more diagnostics and output)
  trackers = update.trackers(stepwise_timer)

  # Code if we want to create a visualization
  if (stepwise_timer > 1000 & stepwise_timer < 1201 & plot.gif == TRUE) {
    ggsave(paste0("Final_Pipeline/Final_Out/ANIM_PLOTS/",
                  stepwise_timer - 1000, ".png"),
           plot = plot.hex.net.gg(),
           width = 6, height = 6, units = "in", dpi = 400)
  }

  # This script allows us to report more diagnostics while the code is running,
  # plot stuff, or even break the loop in a safe way
  source("BreakScript.R")

  # record some other stuff after the biofilm has completely stopped growing
  if (stepwise_timer == 1501) {
    whosig_internal = sample(find.internal(), 500)
    whosig_external = sample(find.external(), 1000)
    whosig_ids = c(whosig_internal, whosig_external)
    whosig_mat = matrix(nrow = ticks - stepwise_timer + 1, ncol = 1500)
  }

  if (stepwise_timer > 1500) {
    whosig_mat[stepwise_timer - 1500, ] = atr[whosig_ids, 7] < atr[whosig_ids, 2]
  }
}

if (track.min.glu == TRUE) {
  min_glu = rbind(occupied, min_glu)
}





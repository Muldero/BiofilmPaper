
### Utility ####################################################################

# Find the most frequently occuring value in a vector
most.freq.val = function(x, k = 1){
  x %>% table() %>% sort(decreasing = TRUE) %>% names() %>% head(k)
}

rescale = function(a, b) { # rescale a to match range of b
  if (typeof(b) == "list") {b = unlist(b)}
  min(b) + (a - min(a)) * (max(b) - min(b)) / (max(a) - min(a))
}

# draw a nonzero value from a row of a vector
draw_nonzero <- function(row) {
  nonzero_values <- row[row != 0]
  if (length(nonzero_values) > 0) {
    return(sample(nonzero_values, 1))
  } else {
    return(0)
  }
}

# Find inherited hunger threshold
draw.hunger.threshold.stdnorm = function(g, herit = thresh_herit) {
  if (sd2 == 0) {
    return(cbind(g, g))
  } else if (herit == 1) {
    g2 = rtnorm(length(g), g, sd2, lower_thresh, upper_thresh)
    return(cbind(g2, g2))
  } else {
    g2 = rtnorm(length(g), g, sd2, lower_thresh, upper_thresh) # new genotype

    # mean of upper and lower
    mean_value = (upper_thresh + lower_thresh)/2
    sd_value = (upper_thresh - lower_thresh)*sqrt(1/12) # sd of uniform
    scaled_g2 = (g2 - mean_value)/sd_value
    inherited = scaled_g2*sqrt(herit)
    deviation = rnorm(length(g), 0, sqrt(1 - herit))
    new_threshhold = (inherited + deviation)*sd_value + mean_value
    return(cbind(new_threshhold, g2))
  }
}

### Growing the biofilm ########################################################

# Create a perfectly hexagonal attribute table (atr)
# Doesn't really fill it, just initializes it to the correct size for our model
make.hex.atr = function(radius = r) { # radius
  atr = matrix(0, nrow = 1 + 3*(radius + 1)*radius, ncol = 17)
  colnames(atr) = c("ID [,1]", "Stress_Threshold [,2]", "ST_Genotype [,3]",
                    "K_e [,4]", "K_i [,5]", "G_e [,6]", "G_i [,7]",
                    "V [,8]", "Depth [,9]", "X [,10]", "Y [,11]",
                    "neigh1 [,12]", "neigh2 [,13]", "neigh3 [,14]",
                    "neigh4 [,15]", "neigh5 [,16]", "neigh6 [,17]")

  atr[, 1] = 1:length(atr[, 1])

  # initialize coordinates and neighbors for first layer
  atr[1, 9:17] = c(0, 0, 0, 2:7)
  atr[2, 9:15] = c(1, 0, 1, 1, 3, 7, NA)
  atr[3, 9:15] = c(1, 1, 0.5, 1, 2, 4, NA)
  atr[4, 9:15] = c(1, 1, -0.5, 1, 3, 5, NA)
  atr[5, 9:15] = c(1, 0, -1, 1, 4, 6, NA)
  atr[6, 9:15] = c(1, -1, -0.5, 1, 5, 7, NA)
  atr[7, 9:15] = c(1, -1, 0.5, 1, 2, 6, NA)

  for (depth in 2:radius) {
    atr = add.layer(depth, at_tab = atr)
  }

  # add in placeholder values for cell parameters
  atr[, 4] = K_m # K_e
  atr[, 5] = 300 # K_i
  atr[, 7] = sample(seq(upper_thresh, G_m, by = 0.01), # G_i
                    length(atr[, 7]), replace = TRUE)
  atr[, 8] = V_0 # V
  if (upper_thresh == lower_thresh) {
    atr[, 2] = upper_thresh # Stress_Thresh
  } else {
    atr[, 2] = sample(seq(lower_thresh, upper_thresh, by = 0.01),
                      length(atr[, 2]), replace = TRUE)
  }
  atr[, 3] = atr[, 2]

  atr[, 12:17][atr[, 12:17] == 0] = NA

  return(atr)
}

# Add a layer to the atr matrix. Called in the previous function
# a lot of math here to make sure neighbors and sutff are correct for a hexagon
add.layer = function(depth, at_tab = atr) {
  if (depth < 2) {
    warning("Depth must be > 1 to add an edge safely")
  }
  # add an edge starting at the point (init_id) and ending one cell before the
  # next point. Depth is the distance from the center (center depth = 0)
  # therefore length of new edge equals depth

  # hex number is the area of the hexagon before the new layer is added
  hex_num = 1 + 3*depth*(depth - 1)
  perimeter_prev = 6*(depth - 1)
  perimeter = 6*depth

  new_ids = hex_num + 1:perimeter

  if (length(at_tab[, 1]) == hex_num) {
    at_tab = rbind(at_tab, matrix(0, ncol = 17, nrow = perimeter))
    at_tab[new_ids, 1] = new_ids
  }

  neighbor_mat = matrix(nrow = perimeter, ncol = 4)

  # go around and calculate neighbors one edge at a time
  for (init_id in ((hex_num + 1) + depth*(0:5))) {
    new_cell_ids = init_id:(init_id + depth - 1)

    edge_num = ceiling((init_id - hex_num)/depth) # (1 = upper right, clockwise)

    if (edge_num == 1) {
      neighbor_mat[init_id - hex_num, ] = c(init_id - perimeter_prev,
                                            init_id + 1,
                                            init_id + perimeter - 1, NA)

      for (i in new_cell_ids[2:depth]) {
        neighbor_mat[i - hex_num, ] = c(i - perimeter_prev + -1:0, i - 1, i + 1)
      }
    } else if (edge_num != 6) {
      neighbor_mat[init_id - hex_num, ] =
        c(init_id - perimeter_prev - edge_num + 1, init_id - 1, init_id + 1, NA)

      for (i in new_cell_ids[2:depth]) {
        neighbor_mat[i - hex_num, ] =
          c(i - perimeter_prev - edge_num + 1 + -1:0, i - 1, i + 1)
      }
    } else {
      neighbor_mat[init_id - hex_num, ] = c(init_id - perimeter_prev - 5,
                                            init_id - 1, init_id + 1, NA)

      if (depth > 2) {
        for (i in new_cell_ids[2:(depth - 1)]) {
          neighbor_mat[i - hex_num, ] = c(i - perimeter_prev - 5 + -1:0,
                                          i - 1, i + 1)
        }
      }

      last_id = new_cell_ids[depth]
      neighbor_mat[last_id - hex_num, ] =
        c(last_id - perimeter - perimeter_prev + 1,
          last_id - perimeter + 0:1, last_id - 1)
    }
  }

  at_tab[new_ids, 12:15] = neighbor_mat

  # update previous layer neighbors
  for (id in (hex_num - perimeter_prev + 1):hex_num) { # prev ids
    if (is.na(at_tab[id, 15])) {
      at_tab[id, 15:17] =
        new_ids[which(neighbor_mat == id, arr.ind = TRUE)[, 1]]
    } else {
      at_tab[id, 16:17] =
        new_ids[which(neighbor_mat == id, arr.ind = TRUE)[, 1]]
    }
  }

  # calculate coordinates for new layer
  x_coords = c(0:(depth - 1), # edge 1
               rep(depth, depth), # edge 2
               depth:1, # edge 3
               0:(1 - depth), # edge 4
               rep(-depth, depth), # edge 5
               -depth:-1) # edge 6

  y_coords = c(depth - (0:(depth - 1))/2, # edge 1
               depth/2 - 0:(depth - 1), # edge 2
               -depth/2 - (0:(depth - 1))/2, # edge 3
               -depth + (0:(depth - 1))/2, # edge 4
               -depth/2 + 0:(depth - 1), # edge 5
               depth/2 + (0:(depth - 1))/2) # edge 6

  at_tab[new_ids, 10] = x_coords
  at_tab[new_ids, 11] = y_coords

  at_tab[new_ids, 9] = depth

  # fill in rest of values
  prev_ids = which(at_tab[, 9] == depth - 1)

  prev_neighbor_mat = neighbor_mat*is.element(neighbor_mat, prev_ids)
  parent_ids = apply(prev_neighbor_mat, 1, draw_nonzero)

  parent_matrix <<- rbind(parent_matrix, cbind(new_ids, parent_ids))

  at_tab[new_ids, 2:8] = at_tab[sample(prev_ids, perimeter, replace = 1), 2:8]
  at_tab[new_ids, 2] = sample(seq(lower_thresh, upper_thresh, by = 0.01),
                              length(new_ids), replace = TRUE)

  at_tab[, 12:17][at_tab[, 12:17] == 0] = NA

  return(at_tab)
}

# Grow a biofilm genetically
# this is what actually fills in the biofilm
grow.a.biofilm2 = function(method = c("entire", "first_layer", "add_layer")) {
  # first we initilize the core of the biofilm
  # just 1 cell occupied
  if (method == "entire" | method == "first_layer") {
    atr = make.hex.atr() # make the matrix
    atr[2:length(atr[, 1]), c(2,3)] = NA
    atr[1, 2:3] = mean(c(upper_thresh, lower_thresh))
    occupied = 1
    border_cells = matrix(c(1, 6), nrow = 1, ncol = 2)
    colnames(border_cells) = c("ID", "Num_Empty_Neighbors")
    parent_matrix = matrix(c(1, 1), nrow = 1, ncol = 2)
    colnames(parent_matrix) = c("Offspring", "Parent")
  }

  # now, if we want to grow the entire thing
  if (method == "entire") {
    while (length(occupied)/length(atr[, 1]) < prop_occupied) {
      # select border cells to reproduce, weight by distance from center
      parent = sample(border_cells[, 1], 1,
                      prob = 1/(atr[border_cells[, 1], 9] + 1))
      potential_daughter = atr[parent, 12:17][
        which(is.element(atr[parent, 12:17], occupied) == FALSE &
                is.na(atr[parent, 12:17]) == FALSE)]
      daughter = potential_daughter[sample(length(potential_daughter), 1)]

      parent_matrix = rbind(parent_matrix, c(daughter, parent))

      atr[daughter, c(2, 3)] = draw.hunger.threshold.stdnorm(atr[parent, 3])
      occupied = c(occupied, daughter)

      border_cells = rbind(border_cells, c(daughter, length(which(is.element(
        atr[daughter, 12:17], occupied) == FALSE &
          is.na(atr[daughter, 12:17]) == FALSE))))

      daughter_occ_neighs = atr[daughter, 12:17][which(is.element(
        atr[daughter, 12:17], occupied))]

      border_cells[which(is.element(border_cells[, 1],
                                    daughter_occ_neighs)), 2] %-=% 1

      border_cells = subset(border_cells, border_cells[, 2] > 0)

      if (length(unique(occupied)) != length(occupied)) {break()}
    }

  # otherwise, if we are just adding a layer
  } else if (method == "add_layer") {
    for (new_babies in 1:ceiling(length(border_cells[, 1])*gens_per_tick)) {
      # choose parents form border cells, weight by distance from center
      parent = sample(border_cells[, 1], 1,
                      prob = 1/(atr[border_cells[, 1], 9] + 1))
      potential_daughter = atr[parent, 12:17][
        which(is.element(atr[parent, 12:17], occupied) == FALSE &
                is.na(atr[parent, 12:17]) == FALSE)]
      daughter = potential_daughter[sample(length(potential_daughter), 1)]

      parent_matrix = rbind(parent_matrix, c(daughter, parent))

      atr[daughter, c(2, 3, 5, 7, 8)] = c(
        draw.hunger.threshold.stdnorm(atr[parent, 3]),
        300, sample(seq(upper_thresh, G_m, by = 0.01), 1), V_0)
      occupied = c(occupied, daughter)

      border_cells = rbind(border_cells, c(daughter, length(which(is.element(
        atr[daughter, 12:17], occupied) == FALSE &
          is.na(atr[daughter, 12:17]) == FALSE))))

      daughter_occ_neighs = atr[daughter, 12:17][which(is.element(
        atr[daughter, 12:17], occupied))]

      border_cells[which(is.element(border_cells[, 1],
                                    daughter_occ_neighs)), 2] %-=% 1

      border_cells = subset(border_cells, border_cells[, 2] > 0)

      if (length(unique(occupied)) != length(occupied)) {break()}
    }
  }
  atr[which(is.element(atr[, 1], occupied) == FALSE), c(5, 7, 8)] = 0
  atr[which(is.element(atr[, 1], occupied) == FALSE), 2:3] = -1

  # this returns to output of this function to the global environment
  # usually not a great idea, but in this case it makes things easy
  parent_matrix <<- parent_matrix
  border_cells <<- border_cells
  occupied <<- occupied
  return(atr)
}

# This function transforms the membrane potential V into the effect on
# glutamate uptake of V. Equation S2 in the paper
Eff = function(V) {
  1/(1 + exp(V - V_0))
}

# Determine signalers
# Larkin et al. (2018) defined signalers after they observed a bimodal
# distribution for membrane potential. My membrane potential is log distributed
# and not bimodal, but F(V) is strongly bimodal. The following function
# identifies signalers by finding the two modes and identifying the cells in
# each one.
find.sigs.fv = function(non = FALSE, external = FALSE, at_tab = atr) {
  if (diff(range(at_tab[occupied, 8])) < 20) {
    # No signalers possible if the biofilm is really small
    return(c())
  } else {
    # The higher mode is always at 1. we want to find the lower mode
    # First, identify a number below this mode and ignore all values
    # greater than it
    F_v = Eff(at_tab[occupied, 8])
    mid_val = weighted.mean(c(mean(F_v), 1), c(2/3, 1/3))

    # Now identify mode for cells below the mode at 1
    lower_half_fv = F_v[which(F_v < mid_val)]
    dens = density(lower_half_fv, adjust = 0.2)
    xval = dens$x[which.max(dens$y)]

    # Now define our cutoff to be a weighted average between this mode and 1
    cutoff = weighted.mean(c(xval, 1), c(3/4, 1/4))

    if (non == FALSE) {
      out = occupied[which(F_v > cutoff)]
    } else { # return nonsignalers
      out = occupied[which(F_v <= cutoff)]
    }

    if (external == FALSE) {
      return(out)
    } else {
      return(out[which(is.element(out, find.external(at_tab)))])
    }
  }
}

# Find cells which are (threshold-defined) signalers
# ie. cells where their internal glutamate is lower than their threshold
find.sigs = function(at_tab = atr) {
  at_tab[, 7] < at_tab[, 2]
}

### Simulation equations #######################################################

# Update potassium across the biofilm - equation S4 in the paper
pot.update = function(at_tab = atr) {

  # Step 1: Signaling (or not)
  sigs = which(find.sigs(at_tab) == 1)
  signal_kplus = F_param*(g_k/(tpp/2))*(at_tab[sigs, 8] -
                                      100*log(at_tab[sigs, 4]/at_tab[sigs, 5]))

  at_tab[sigs, 4] = at_tab[sigs, 4] + signal_kplus*(signal_kplus > 0)
  at_tab[sigs, 5] = at_tab[sigs, 5] - signal_kplus*(signal_kplus > 0)

  # Step 2: Uptake
  uptake_kplus = (D_p/(tpp/2))*at_tab[, 4]*at_tab[, 7]*(300 - at_tab[, 5])
  at_tab[, 4] = at_tab[, 4] - uptake_kplus
  at_tab[, 5] = at_tab[, 5] + uptake_kplus

  # if k_e negative then need to subtract from k_i. There's no way for k+ to
  # enter the system, so need it to not be lost either...
  negs = which(at_tab[, 4] < 0)
  if (length(negs > 0)) {
    at_tab[negs, 5] = at_tab[negs, 5] + at_tab[negs, 4]
    at_tab[negs, 4] = 0
  }

  # average potassium among occupied cells to simulate diffusion
  at_tab[occupied, 4] = mean(c(at_tab[occupied, 4]))

  return(at_tab)
}

# Update glutamate across the biofilm - equation S3
glu.update = function(at_tab = atr) {

  # Diffusion
    # begin at the edges
  edges = which(at_tab[, 9] == max(at_tab[, 9])) # edgeid
  rank_ids = edges

  orig_glu = sum(at_tab[rank_ids, 7])

  # glutamate uptake by cells
  at_tab[rank_ids, 6:7] = glu.calc.uptake(G_m_vec[1:length(rank_ids)],
                                               rank_ids, at_tab)

  # this will let us keep track of unabsorbed glu after everything has filtered
  absorbed_glu = sum(at_tab[rank_ids, 7]) - orig_glu

  # moving in from the outside to the center of the biofilm
  for (rank in (max(at_tab[, 9]) - 1):0) {
    prev_ids = rank_ids # cells one rank out
    rank_ids = which(at_tab[, 9] == rank)

    # diffusion
    neighbors = at_tab[rank_ids, 12:17]
    neighbors = matrix(ifelse(is.element(neighbors, prev_ids),
                              at_tab[neighbors, 6], NA), ncol = 6)
    avail_glu = rowMeans(neighbors, na.rm = TRUE)

    orig_glu = at_tab[rank_ids, 7] # glu before uptake
    # uptake by cells
    at_tab[rank_ids, 6:7] = glu.calc.uptake(avail_glu, rank_ids, at_tab)
    absorbed_glu = absorbed_glu + sum(at_tab[rank_ids, 7] - orig_glu) # after
  }

  # This is letting glu carry over to the next iteration
  at_tab[, 6] = max((G_m*length(edges) - absorbed_glu)/length(at_tab[, 1]), 0)

  # metabolism
  at_tab[, 7] = at_tab[, 7] - (delta_g/(tpp/2))*at_tab[, 7]

  # return output
  return(at_tab)
}

# Calculate glutamate uptake for an individual cell equation S1
glu.calc.uptake = function(avail_vec, ids, at_tab) {
  uptake_vec = (alpha_g/(tpp/2))*Eff(at_tab[ids, 8])*
    avail_vec/(0.75 + avail_vec)

  # Check to see if can actually uptake as much as possible or not
  unlimited = avail_vec > uptake_vec
  if (length(unlimited) > 0) {
    avail_vec[unlimited] = avail_vec[unlimited] -
      uptake_vec[unlimited]
    at_tab[ids[unlimited], 7] = at_tab[ids[unlimited], 7] +
      uptake_vec[unlimited]
  }
  limited = avail_vec <= uptake_vec
  if (length(limited) > 0) {
    at_tab[ids[limited], 7] = at_tab[ids[limited], 7] + avail_vec[limited]
    avail_vec[limited] = 0
  }

  # update available glu for next rank of cells
  at_tab[ids, 6] = avail_vec

  # return important numbers
  return(at_tab[ids, 6:7])
}

# Update membrane potential V - equation S5
mem.potential.update = function(at_tab = atr) {
  leak = -156 + 4*(at_tab[, 4] - K_m)/(1 - exp(10*(K_m - at_tab[, 4])))
  leak[which(is.na(leak))] = -156 + 4*0.1

  Nernst_pot = 100 * log(at_tab[, 4]/at_tab[, 5])
  Nernst_pot[which(Nernst_pot == Inf)] = 0

  delta = -(18/(tpp/2))*(at_tab[, 8] - leak) -
    (g_k/(tpp/2))*(at_tab[, 8] - Nernst_pot)*find.sigs(at_tab)*
    ((g_k/(tpp/2))*(at_tab[, 8] - Nernst_pot) > 0)

  at_tab[, 8] = at_tab[, 8] + delta

  return(at_tab)
}

### Plotting ###################################################################

potfunc = colorRampPalette(c("palegreen", "black"))
pot_grad = potfunc(160)
glufunc = colorRampPalette(c("lightpink", "black"))
glu_grad = glufunc(490)
thrfunc = colorRampPalette(c("lightyellow", "black"))
thr_grad = thrfunc(360)
voltfunc = colorRampPalette(c("black", "cyan"))
volt_grad = voltfunc(101)

# Plot the biofilm
plot.hex.net.gg = function(remake = FALSE, # recalculate plotting coordinates?
                       return.color.df = FALSE, # don't plot, just give colors
                       at_tab = atr,
                       cat = NA, altcols = NA, # can deifne alternate colors
                       are.unoccupied = NA, # define unoccupied explicitely
                       plot.boundary = FALSE # plot inner/outer boundary
                       ) {

  # This weird scaling reflects the actual effect of V on K+ uptake
  voltage_colors = floor(abs(at_tab[, 8])*10 - 1519)
  voltage_colors[which(voltage_colors > 81)] = 81
  voltage_colors[which(voltage_colors < 1)] = 1
  frame_color1 = (volt_grad[c(21,101)])[voltage_colors]
  point_color1 = (volt_grad[c(1:81)])[voltage_colors]

  # use alternate colors?
  if (length(altcols) > 1) {
      colors = altcols[cat + max(cat) + 1]

    # explicitely defined unoccupied cells?
    if (length(are.unoccupied) > 1) {
      colors[are.unoccupied] = "#FFFFFF"
      occupied = which(1:length(at_tab[, 1]) %in% are.unoccupied == FALSE)
    }

    frame_color1 = colors
    point_color1 = colors
  }

  # plot inner vs outer boundary?
  if (plot.boundary == TRUE) {
    bounds = which(atr[, 9] %in% (find.glu.radius() + c(-1, 0, 1)))
    frame_color1[bounds] = "yellow"
    point_color1[bounds] = "yellow"
  }

  # calculate plotting coordinates for cells
  # Will automatically do if code ahsn't been run before, otherwise needs to be
  # told to recalculate
  if (remake == "TRUE" | exists("hexdat") == FALSE) {
    lx = at_tab[, 10]
    ly = at_tab[, 11]
    hex_x = cbind(
      lx - 1/(2*sqrt(3)) - 0.06699,
      lx - 1/(sqrt(3)) - 0.06699,
      lx - 1/(2*sqrt(3)) - 0.06699,
      lx + 1/(2*sqrt(3)) + 0.06699,
      lx + 1/(sqrt(3)) + 0.06699,
      lx + 1/(2*sqrt(3)) + 0.06699
    )
    hex_y = cbind(ly + 0.5, ly + 0, ly - 0.5, ly - 0.5, ly + 0, ly + 0.5)
    hex_vectors <<- cbind(hex_x, hex_y)

    # Get coordinates in long format with an id
    hexdat_x <- reshape2::melt(cbind(id = 1:length(at_tab[, 1]),
                                     as.data.frame(hex_x)),
                     id.vars = "id", value.name = "x")
    hexdat_y <- reshape2::melt(cbind(id = 1:length(at_tab[, 1]),
                                     as.data.frame(hex_y)),
                     id.vars = "id", value.name = "y")

    # Merge them into the same dataframe
    hexdat <- merge(hexdat_x, hexdat_y)
  }

  colors = cbind(id = at_tab[, 1], col = point_color1)
  plot_df = merge(hexdat, colors)
  plot_df = plot_df[plot_df$id %in% occupied, ]

  # just return to colors dataframe instead of plotting it?
  if (return.color.df) {
    return(plot_df)
  }

  # plot it!
  p1 = ggplot(plot_df, aes(x, y)) +
    geom_polygon(aes(group = id), fill = plot_df$col, colour = "#00000000") +
    coord_fixed(ratio = max(at_tab[occupied, 10])/max(at_tab[occupied, 11])) +
    theme_void()

  return(p1)
}

# alternate plotting code in baseR
plot.hex.net.quickndirty = function(remake = FALSE,
                                    as.hex = TRUE, at_tab = atr,
                                    plot.boundary = FALSE) {

  # This weird scaling reflects the actual effect of V on K+ uptake
  voltage_colors = floor(abs(at_tab[, 8])*10 - 1519)
  voltage_colors[which(voltage_colors > 81)] = 81
  voltage_colors[which(voltage_colors < 1)] = 1
  frame_color1 = (volt_grad[c(21,101)])[voltage_colors]
  point_color1 = (volt_grad[c(1:81)])[voltage_colors]

  if (plot.boundary == TRUE) {
    bounds = which(at_tab[, 9] %in% (find.glu.radius() + c(-1, 0, 1)))
    frame_color1[bounds] = "yellow"
    point_color1[bounds] = "yellow"
  }

  if (remake == "TRUE" | exists("hex_vectors") == FALSE) {
    lx = at_tab[, 10]
    ly = at_tab[, 11]
    hex_x = cbind(
      lx - 1/(2*sqrt(3)) - 0.06699,
      lx - 1/(sqrt(3)) - 0.06699,
      lx - 1/(2*sqrt(3)) - 0.06699,
      lx + 1/(2*sqrt(3)) + 0.06699,
      lx + 1/(sqrt(3)) + 0.06699,
      lx + 1/(2*sqrt(3)) + 0.06699
    )
    hex_y = cbind(ly + 0.5, ly + 0, ly - 0.5, ly - 0.5, ly + 0, ly + 0.5)
    hex_vectors <<- cbind(hex_x, hex_y)
  }

  # Signaling
  point_color1[unoccupied] = "white"

  plot(1, col = "white", xlim = c(-r, r), ylim = c(r, -r),
       xlab = "", ylab = "", axes = FALSE, frame.plot = FALSE, asp = 1,
       main = "", col.main = "white")

  polygon(x = c(0, -r - 0.64, -r - 0.64, 0, r + 0.64, r + 0.64),
          y = c(min(at_tab[, 11]) - 0.8,
                -r/2 - 0.5,
                r/2 + 0.5,
                max(at_tab[, 11]) + 0.8,
                r/2 + 0.5,
                -r/2 - 0.5),
          col = "white", border = NA)


  for (i in occupied){#not_common) {
    polygon(x = hex_vectors[i, 1:6], y = hex_vectors[i, 7:12],
            border = NA, col = point_color1[i], lty = par("lty"),
            fillOddEven = TRUE)
  }
}

# Return a plot of the biofilm as a .tiff file
# much much faster than having it actually show up in R
plot.a.tiff.please = function(title,
                          location = "Plot_Code/Final_Plots/Images/",
                          plot.boundary = FALSE) {
  tiff(paste0(location, title, ".tiff"), width = 5.2, height = 5.2,
       units = "in", res = 600)
  plot.hex.net.quickndirty(plot.boundary = plot.boundary)
  dev.off()
}

### Diagnostic/Analysis Functions ##############################################

# Midline is some y values that every oscillation crosses. This just finds the
# maxima between each location where we cross the midline.
# To use correctly, plot the trace, look at the plot, and manually identify a
# midline value that every single oscillation crosses.
#
# Check results! This code is not 100% perfect!
get.only.max.peaks = function(v, midline = 0.15, top = TRUE) {
  adjusted_v = v - midline

  pos_cross = c()
  max_ids = c()

  for (i in 2:length(adjusted_v)) {
    if (adjusted_v[i - 1] <= 0 & adjusted_v[i] > 0 & top == TRUE) {
      pos_cross = c(pos_cross, i)
      lpc = length(pos_cross)
      if (lpc > 1) {
        max_ids = c(max_ids, pos_cross[lpc - 1] - 1 +
                      which.max(adjusted_v[pos_cross[lpc - 1]:pos_cross[lpc]]))
      }
    }

    if (adjusted_v[i - 1] <= 0 & adjusted_v[i] > 0 & top == FALSE) {
      pos_cross = c(pos_cross, i)
      lpc = length(pos_cross)
      if (lpc > 1) {
        max_ids = c(max_ids, pos_cross[lpc - 1] - 1 +
                      which.min(adjusted_v[pos_cross[lpc - 1]:pos_cross[lpc]]))
      }
    }
  }

  return(max_ids)
}


### Identify internal vs external cells
# Defined as the depth at which glutamate doesn't filter if every cell can
# absorb to maximum capacity
find.glu.radius = function(at_tab = atr) {
  at_tab_test = at_tab
  at_tab_test[occupied, 8] = -156
  at_tab_test[occupied, 4] = K_m
  at_tab_test[occupied, 5] = 300
  old_glu = at_tab_test[occupied, 7]

  ## Now run to determine glu uptaken for each cell
  edges = which(at_tab_test[, 9] == max(at_tab_test[, 9])) # edgeid
  rank_ids = edges
  orig_glu = sum(at_tab_test[rank_ids, 7])
  # external glu = G_m
  at_tab_test[rank_ids, 6:7] = glu.calc.uptake(G_m_vec[1:length(rank_ids)],
                                               #rep(G_m, length(rank_ids)),
                                               rank_ids, at_tab_test)
  absorbed_glu = sum(at_tab_test[rank_ids, 7]) - orig_glu
  for (rank in (max(at_tab_test[, 9]) - 1):0) {
    prev_ids = rank_ids # cells one rank out
    rank_ids = which(at_tab_test[, 9] == rank)
    neighbors = at_tab_test[rank_ids, 12:17]
    neighbors = matrix(ifelse(is.element(neighbors, prev_ids),
                              at_tab_test[neighbors, 6], NA), ncol = 6)
    avail_glu = rowMeans(neighbors, na.rm = TRUE)
    orig_glu = at_tab_test[rank_ids, 7] # glu before uptake
    at_tab_test[rank_ids, 6:7] =
      glu.calc.uptake(avail_glu, rank_ids, at_tab_test)
    absorbed_glu =
      absorbed_glu + sum(at_tab_test[rank_ids, 7] - orig_glu) # after
  }
  at_tab_test[, 6] = max((G_m*length(edges) -
                            absorbed_glu)/length(at_tab_test[, 1]), 0)

  new_glu = at_tab_test[occupied, 7]
  updated = at_tab_test[occupied, ]
  updated[, 7] = updated[, 7] - old_glu

  depth_means = tapply(updated[, 7], updated[, 9], mean)
  # plot(depth_means, ylim = c(0, 1))
  if (sum(depth_means == 0) == 0) {
    return(0)
  } else {
    return(min(which(depth_means > 0)))
  }
}

# determine occupied internal cells
find.internal = function(at_tab = atr) {
  int_indx = at_tab[which(at_tab[, 9] < glu_radius), 1]
  return(int_indx[which(is.element(int_indx, occupied))])
}

# determine occupied external cells
find.external = function(at_tab = atr) {
  # radius = max(at_tab[occupied, 9])
  ext_indx = at_tab[which(at_tab[, 9] >= glu_radius), 1]
  return(ext_indx[which(is.element(ext_indx, occupied))])
}

### Record who is signaling each iteration
check.whos.signaling = function(t_tab = trackers, at_tab = atr) {
  if (is.growing == FALSE) {
    temp_sigs = find.sigs.fv()
    if (length(temp_sigs)/length(occupied) >=
        trackers$Tot_fv_sig[stepwise_timer - 1]) {

      who_is_signaling[, length(who_is_signaling[1, ])] =
        c(occupied %in% temp_sigs)

      colnames(who_is_signaling)[length(colnames(who_is_signaling))] =
        stepwise_timer

    } else if (colnames(who_is_signaling)[length(colnames(who_is_signaling))] ==
               (stepwise_timer - 1)) {
      who_is_signaling = cbind(who_is_signaling, c(occupied %in% temp_sigs))
      colnames(who_is_signaling)[length(colnames(who_is_signaling))] = 0
    }
  }
  return(who_is_signaling)
}


### Record every parameter we are actually interested in for final analysis

# In: inner, Out: outer, Tot: total, Sig: signalers (defined by fv),
# Non: nonsignalers, Thr: threshold, Inher: inheritance rate,
# Var: variance, Med: median, all else mean
trackers = matrix(0, ncol = 38, nrow = ticks)
colnames(trackers) = c("Time", "Radius",
                       "In_fv_sig", "In_sig", "In_sig_thr", "In_non_thr",
                       "In_non_fv", "In_med_V", "In_G_i", "In_K_e", "In_K_i",
                       "In_var_G_i", "In_sig_inher", "In_non_inher",
                       "Out_fv_sig", "Out_sig", "Out_sig_thr", "Out_non_thr",
                       "Out_non_fv", "Out_med_V", "Out_G_i", "Out_K_e",
                       "Out_K_i", "Out_var_G_i", "Out_sig_inher",
                       "Out_non_inher",
                       "Tot_fv_sig", "Tot_sig", "Tot_sig_thr", "Tot_non_thr",
                       "Tot_non_fv", "Tot_med_V", "Tot_G_i", "Tot_K_e",
                       "Tot_K_i", "Tot_var_G_i", "Tot_sig_inher",
                       "Tot_non_inher")

trackers = trackers %>% as.data.frame()

update.trackers = function(iteration, t_tab = trackers, at_tab = atr) {
  # extend trackers if too short
  if (length(t_tab[, 1]) < iteration) {
    t_tab = rbind(t_tab, t_tab[1, ])
  }

  t_tab[iteration, ] = c(iteration, max(at_tab[occupied, 9]),
                         get.measure.tracker(find.internal(at_tab), at_tab),
                         get.measure.tracker(find.external(at_tab), at_tab),
                         get.measure.tracker(occupied, at_tab)
  )
  return(t_tab)
}

get.measure.tracker = function(ids, at_tab = atr) {
  signalers = ids[which(is.element(ids, find.sigs.fv()))]
  nonsignalers = ids[which(is.element(ids, signalers) == FALSE)]

  parent_matrix_t = parent_matrix[which(is.element(parent_matrix[, 1],
                                                   ids)), ]
  parent_matrix_t = parent_matrix_t[which(is.element(parent_matrix_t[, 2],
                                                     ids)), ]
  daughter_states = is.element(parent_matrix_t[, 1], signalers)
  parent_states = is.element(parent_matrix_t[, 2], signalers)
  sig_inher_rate = sum(daughter_states == parent_states &
                         daughter_states == 1)/sum(daughter_states == 1)
  non_inher_rate = sum(daughter_states == parent_states &
                         daughter_states == 0)/sum(daughter_states == 0)

  return(c(
    length(signalers)/length(ids),                         # fv_sig
    sum(at_tab[ids, 7] <= at_tab[ids, 2])/length(ids),     # sig
    mean(at_tab[signalers, 2]),                            # sig_thr
    mean(at_tab[nonsignalers, 2]),                         # non_thr
    mean(Eff(at_tab[nonsignalers, 8])),                    # non_fv
    median(at_tab[ids, 8]),                                # med_V
    mean(at_tab[ids, 7]),                                  # G_i
    mean(at_tab[ids, 4]),                                  # K_e
    mean(at_tab[ids, 5]),                                  # K_i
    var(at_tab[ids, 7]),                                   # var_G_i
    sig_inher_rate,                                        # sig_inher
    non_inher_rate                                         # non_inher
  ))
}

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

# Create a perfectly hexagonal attribute table (atr)
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

# Add a layer to the atr
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

  if (do.inheritance == TRUE) {
    at_tab[new_ids, 4:8] = at_tab[parent_ids, 4:8]
    at_tab[new_ids, 2:3] = draw.hunger.threshold.stdnorm(at_tab[parent_ids, 3])
  } else {
    at_tab[new_ids, 2:8] = at_tab[sample(prev_ids, perimeter, replace = 1), 2:8]
    at_tab[new_ids, 2] = sample(seq(lower_thresh, upper_thresh, by = 0.01),
                                length(new_ids), replace = TRUE)
  }

  at_tab[, 12:17][at_tab[, 12:17] == 0] = NA

  return(at_tab)
}

# Grow a biofilm genetically
grow.a.biofilm2 = function(method = c("entire", "first_layer", "add_layer")) {
  if (method == "entire" | method == "first_layer") {
    atr = make.hex.atr()
    atr[2:length(atr[, 1]), c(2,3)] = NA
    atr[1, 2:3] = mean(c(upper_thresh, lower_thresh))
    occupied = 1
    border_cells = matrix(c(1, 6), nrow = 1, ncol = 2)
    colnames(border_cells) = c("ID", "Num_Empty_Neighbors")
    parent_matrix = matrix(c(1, 1), nrow = 1, ncol = 2)
    colnames(parent_matrix) = c("Offspring", "Parent")
  }

  if (method == "entire") {
    while (length(occupied)/length(atr[, 1]) < prop_occupied) {
      parent = sample(border_cells[, 1], 1,
                      prob = 1/(atr[border_cells[, 1], 9] + 1))
      # 1/border_cells[, 2]) # weight by neighbor
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

      if (length(occupied) %% 1000 == 1) {
        # flexy.plot(is.element(atr[, 1], occupied), at_tab = atr)
        print(length(occupied))
      }
    }
  } else if (method == "add_layer") {
    for (new_babies in 1:ceiling(length(border_cells[, 1])*gens_per_tick)) {
      # limit how spread out borders can be
      # elig_border_cells = border_cells[which(atr[border_cells[, 1], 9] <=
      #                     min(atr[border_cells[, 1], 9]) + 5), , drop = FALSE]
      # parent = sample(elig_border_cells[, 1], 1,
      #                 prob = 1/elig_border_cells[, 2]) # weight by neighbor
      parent = sample(border_cells[, 1], 1,
                      prob = 1/(atr[border_cells[, 1], 9] + 1))
      # prob = 1/border_cells[, 2]) # weight by neighbor
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

      if (length(occupied) %% 1000 == 1) {
        # flexy.plot(is.element(atr[, 1], occupied), at_tab = atr)
        # flexy.plot(floor((atr[occupied, 2])*100), at_tab = atr[occupied,])
      }
      # source("BreakScript.R")
    }
  }
  atr[which(is.element(atr[, 1], occupied) == FALSE), c(5, 7, 8)] = 0
  atr[which(is.element(atr[, 1], occupied) == FALSE), 2:3] = -1

  parent_matrix <<- parent_matrix
  border_cells <<- border_cells
  occupied <<- occupied
  return(atr)
}


Eff = function(V) {
  1/(1 + exp(V - V_0))
}

calc.benefit.atr = function(type = c("individual", "group"), at_tab = atr) {

  nonsig_fv = mean(Eff(at_tab[find.sigs.fv(non = TRUE), 8]))
  sig_fv = mean(Eff(at_tab[find.sigs.fv(non = FALSE), 8]))

  if (type == "group") {
    ns_avail_glu = at_tab[find.sigs.fv(non = TRUE), 6]
    s_avail_glu = at_tab[find.sigs.fv(non = FALSE), 6]

    uptake_ns = (alpha_g/(tpp/2))*nonsig_fv*ns_avail_glu/(0.75 + ns_avail_glu)
    uptake_s = (alpha_g/(tpp/2))*sig_fv*s_avail_glu/(0.75 + s_avail_glu)

    return(c(mean(uptake_s), mean(uptake_ns)))
  } else {
    # col1 = sig, col2 = nonsig
    out = matrix(0, nrow = length(atr[, 1]), ncol = 2)

    fv = Eff(at_tab[, 8])
    avail_glu = at_tab[, 6]

    uptake = (alpha_g/(tpp/2))*fv*avail_glu/(0.75 + avail_glu)
    alt_uptake = (alpha_g/(tpp/2))*nonsig_fv*avail_glu/(0.75 + avail_glu)*
                    is.element(at_tab[, 1], find.sigs.fv()) +
                 (alpha_g/(tpp/2))*sig_fv*avail_glu/(0.75 + avail_glu)*
                    is.element(at_tab[, 1], find.sigs.fv(non = TRUE))

    ns_id = find.sigs.fv(non = TRUE)
    s_id = find.sigs.fv()

    out[s_id, 1] = uptake[s_id]
    out[s_id, 2] = alt_uptake[s_id]
    out[ns_id, 2] = uptake[ns_id]
    out[ns_id, 1] = alt_uptake[ns_id]

    return(out)
  }
}

# Determine signalers
# Joe's paper defined signalers after they observed a bimodal distribution for
# membrane potential. My membrane potential is log distributed and not bimodal,
# F(V) is strongly bimodal. The following function identifies signalers by
# finding the two modes and identifying the cells in each one.
#
# Note that this works most of the time but every so often it gets confused.
# However, when that happens it still gets something in between the two maxima.
find.sigs.fv = function(non = FALSE, external = FALSE, at_tab = atr) {
  if (diff(range(at_tab[occupied, 8])) < 20) {
    return(c())
  } else {
    F_v = Eff(at_tab[occupied, 8])
    mid_val = weighted.mean(c(mean(F_v), 1), c(2/3, 1/3))

    # Now identify the maxima for cells below the midpoint
    lower_half_fv = F_v[which(F_v < mid_val)]
    dens = density(lower_half_fv, adjust = 0.2)
    xval = dens$x[which.max(dens$y)]

    # Now define our cutoff to be a weighted average between x and midval
    cutoff = weighted.mean(c(xval, 1), c(3/4, 1/4))

    # }
    if (non == FALSE) {
      out = occupied[which(F_v > cutoff)]
    } else { # return nonsignalers
      out = occupied[which(F_v <= cutoff)]
    }
    # } else { # no clear distinction
    #   if (non == FALSE) {
    #     out = occupied[which(at_tab[occupied, 7] < at_tab[occupied, 2])]
    #   } else { # return nonsignalers
    #     out = occupied[which(at_tab[occupied, 7] >= at_tab[occupied, 2])]
    #   }
    # }
    # if (length(out)/length(occupied) > 0.95 & non == FALSE) {
    #   mid_val = (median(F_v) + max(F_v))/2
    #   out = occupied[which(F_v > mid_val)]
    # }
    if (external == FALSE) {
      return(out)
    } else {
      return(out[which(is.element(out, find.external(at_tab)))])
    }
  }
}

# Find cells which are (threshold-defined) signalers and not refractory
find.sigs = function(at_tab = atr) {
  at_tab[, 7] < at_tab[, 2]
}

# Calculate fitness of a cell
# Above a certain threshold of G_i the fitness should be 1, at very low G_i it
# should be near 0
get.fit = function(ids, at_tab = atr) {
  gi = at_tab[ids, 7]

  # calculate fitness assuming linear scaling between 1 and min_growth
  # for max_gi to min_gi
  return(((gi - min_gi)/max_gi + 0.4)/1.4)
}

# Check cluster sizes in a biofilm
check.clustersize = function(ids, at_tab = atr) {
  signalers = ids[which(is.element(ids, find.sigs.fv()))]
  sizes = c()
  while (length(signalers > 0)) {
    cluster_ingroupid = signalers[1]
    justadds = cluster_ingroupid
    while (length(justadds) != 0) {
      neighs = c()
      for (i in justadds) {
        neighs = sort(unique(c(neighs, at_tab[i, 12:17])))
      }
      neighs = subset(neighs, neighs > 0)
      neighs = subset(neighs, is.element(neighs, signalers))
      neighs = subset(neighs, is.element(neighs, cluster_ingroupid) == FALSE)
      justadds = neighs
      cluster_ingroupid = c(cluster_ingroupid, justadds)
    }
    sizes = sort(c(sizes, length(cluster_ingroupid)))
    signalers = subset(signalers, is.element(signalers, cluster_ingroupid) == FALSE)
  }
  return(sizes)
}
### Concentration calculations #################################################

# Update potassium across the biofilm
pot.update = function(at_tab = atr) {

  # Step 1: Signaling?
  sigs = which(find.sigs(at_tab) == 1)
  signal_kplus = F_param*(g_k/(tpp/2))*(at_tab[sigs, 8] -
                                      100*log(at_tab[sigs, 4]/at_tab[sigs, 5]))
  # signal_kplus[which(signal_kplus < 0)] = 0
  # print(length(which(signal_kplus < 0))/length(signal_kplus))
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

  # warning("Experimental potassium averaging code in use")
  # if (is.growing == TRUE) {
    at_tab[occupied, 4] = mean(c(at_tab[occupied, 4]))#,
                                 # rep(K_m, floor(length(occupied)*kflux))))
    # rep(K_m, 6*max(at_tab[occupied, 9]))))
  # } else {
  #   at_tab[, 4] = mean(at_tab[, 4])
  # }

  return(at_tab)
}

# Update glutamate across the biofilm
glu.update = function(at_tab = atr) {
  edges = which(at_tab[, 9] == max(at_tab[, 9])) # edgeid
  rank_ids = edges

  orig_glu = sum(at_tab[rank_ids, 7])
  # external glu = G_m
  at_tab[rank_ids, 6:7] = glu.calc.uptake(G_m_vec[1:length(rank_ids)],
                                               #rep(G_m, length(rank_ids)),
                                               rank_ids, at_tab)
  # this will let us keep track of unabsorbed glu after everything has filtered
  absorbed_glu = sum(at_tab[rank_ids, 7]) - orig_glu

  for (rank in (max(at_tab[, 9]) - 1):0) {
    prev_ids = rank_ids # cells one rank out
    rank_ids = which(at_tab[, 9] == rank)

    neighbors = at_tab[rank_ids, 12:17]
    neighbors = matrix(ifelse(is.element(neighbors, prev_ids),
                              at_tab[neighbors, 6], NA), ncol = 6)
    avail_glu = rowMeans(neighbors, na.rm = TRUE)

    orig_glu = at_tab[rank_ids, 7] # glu before uptake
    at_tab[rank_ids, 6:7] = glu.calc.uptake(avail_glu, rank_ids, at_tab)
    absorbed_glu = absorbed_glu + sum(at_tab[rank_ids, 7] - orig_glu) # after
  }

  # This is letting glu carry over to the next iteration. Hopefully keeps
  # inner cells less stressed
  at_tab[, 6] = max((G_m*length(edges) - absorbed_glu)/length(at_tab[, 1]), 0)

  # degrade glutamate
  at_tab[, 7] = at_tab[, 7] - (delta_g/(tpp/2))*at_tab[, 7]

  # return output
  return(at_tab)
}

# Calculate glutamate uptake for an individual cell
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

# Update membrane potential V
mem.potential.update = function(at_tab = atr) {
  leak = -156 + 4*(at_tab[, 4] - K_m)/(1 - exp(10*(K_m - at_tab[, 4])))
  leak[which(is.na(leak))] = -156 + 4*0.1

  Nernst_pot = 100 * log(at_tab[, 4]/at_tab[, 5])
  Nernst_pot[which(Nernst_pot == Inf)] = 0

  delta = -(18/(tpp/2))*(at_tab[, 8] - leak) -
    (g_k/(tpp/2))*(at_tab[, 8] - Nernst_pot)*find.sigs(at_tab)*
    ((g_k/(tpp/2))*(at_tab[, 8] - Nernst_pot) > 0)

  at_tab[, 8] = at_tab[, 8] + delta

  # print(length(which((70/(tpp/2))*(at_tab[, 8] -
  # Nernst_pot)*find.sigs(at_tab) < 0))/length(find.sigs(at_tab)))

  return(at_tab)
}

### Alternate Glu update method (directional uptake)

glu.update.directional = function(at_tab = atr) {
  store_col = at_tab[, 6:7]*0
  at1 = at_tab

  # for (d in 1:6) {
  #   store_col = store_col + filter_glu_by_dir(at_tab, d)
  # }

  store_col = cbind(filter_glu_by_dir(at_tab, 1)[, 2],
                    filter_glu_by_dir(at_tab, 2)[, 2],
                    filter_glu_by_dir(at_tab, 3)[, 2],
                    filter_glu_by_dir(at_tab, 4)[, 2],
                    filter_glu_by_dir(at_tab, 5)[, 2],
                    filter_glu_by_dir(at_tab, 6)[, 2])

  at_tab = at1

  # at_tab[, 6:7] = at_tab[, 6:7] + store_col

  at_tab[, 7] = at_tab[, 7] + apply(store_col, 1, new.mean)

  # degrade glutamate
  at_tab[, 7] = at_tab[, 7] - (delta_g/(tpp/2))*at_tab[, 7]

  # return output
  return(at_tab)
}

filter_glu_by_dir = function(at_tab, direction) {
  init_vals = at_tab[, 6:7]

  edges = which(depth_mat[, direction] == 1) # edgeid
  rank_ids = edges

  orig_glu = sum(at_tab[rank_ids, 7])
  # external glu = G_m
  at_tab[rank_ids, 6:7] = glu.calc.uptake(
    G_m_vec[r*(direction - 1) + 1:(r + 1)], rank_ids, at_tab)
  # 1:length(rank_ids), rep(G_m, length(rank_ids)),

  # this will let us keep track of unabsorbed glu after everything has filtered
  absorbed_glu = sum(at_tab[rank_ids, 7]) - orig_glu

  # rank_ids_full = rank_ids

  for (rank in 2:max(depth_mat[, 1])) {
    prev_ids = rank_ids # _full # cells one rank out
    rank_ids = which(depth_mat[, direction] == rank)

    neighbors = at_tab[rank_ids, 12:17]
    neighbors = matrix(ifelse(is.element(neighbors, prev_ids),
                              at_tab[neighbors, 6], NA), ncol = 6)
    avail_glu = rowMeans(neighbors, na.rm = TRUE)

    orig_glu = at_tab[rank_ids, 7] # glu before uptake
    at_tab[rank_ids, 6:7] = glu.calc.uptake(avail_glu, rank_ids, at_tab)
    absorbed_glu = absorbed_glu + sum(at_tab[rank_ids, 7] - orig_glu) # after
  }

  # This is letting glu carry over to the next iteration. Hopefully keeps
  # inner cells less stressed
  at_tab[, 6] = max((G_m*length(edges) - absorbed_glu)/length(at_tab[, 1]), 0)

  return(at_tab[, 6:7] - init_vals)
}

### Bumping Functions ##########################################################

# Find midpoint assuming using the hex grid system
find.midpoint.xy = function(coord1, coord2) {
  x1 = coord1[1]
  y1 = coord1[2]
  x2 = coord2[1]
  y2 = coord2[2]

  if (abs(x1 - x2) > 1) {
    midx = sample(c(floor(mean(c(x1, x2))), ceiling(mean(c(x1, x2)))), 1)
  } else {
    midx = sample(c(x1, x2), 1)
  }

  if (abs(y1 - y2) >= 1) {
    midy = round(mean(c(y1, y2)))
    if (midx %% 2 == 0) {
      midy = midy + sample(c(0.5, -0.5), 1)
    }
  } else {
    midy = sample(c(y1, y2), 1)
    if (y1 == y2) {
      midy = y1 + sample(c(0.5, -0.5), 1)
    }
    if (midx %% 2 == 0) { # x even
      if (midy %% 1 != 0.5) { # and y integer
        if (midy == y1) {
          midy = y2
        } else {
          midy = y1
        }
      }
    } else if (midy %% 1 != 0) { # x odd but y not integer
      if (midy == y1) {
        midy = y2
      } else {
        midy = y1
      }
    }
  }

  if (midy == 0) {
    midy = 1
  } else if (midy == M + 0.5) {
    midy = M - 0.5
  }

  return(matrix(c(midx, midy), 1, 2))
}

# turn x into a row matrix
arm = function(x) { # as.row.matrix
  matrix(x, nrow = 1, ncol = length(x))
}

# left_path is a matrix 2 columns each row is a coordinate pair
# same for right path. recursively repeat until bottom cell of left and top of
# right are the same or adjacent
# This produces a slightly noisy path. I love it!
shortest.path.recursive = function(leftpath, rightpath, at_tab = atr) {
  if (abs(leftpath[nrow(leftpath), 1] - rightpath[1, 1]) > 1 |
      abs(leftpath[nrow(leftpath), 2] - rightpath[1, 2]) > 1) {
    midpoint = find.midpoint.xy(leftpath[nrow(leftpath), 1:2],
                                rightpath[1, 1:2])
    leftpath = rbind(leftpath, shortest.path.recursive(leftpath, midpoint))
    rightpath = rbind(shortest.path.recursive(midpoint, rightpath), rightpath)
  }
  output = rbind(leftpath, rightpath)
  return(output[!duplicated(output),])
}

# Given a dead cell and the new parent, find the shortest line between them,
# bump all cells in the line down, and create a daughter cell at the start
# (adjacent to the parent). Then recalculate connection status for all cells
# This one is specific to the hunger model
bump.cells.hunger = function(parentid, deadid, at_tab = atr,
                             sh_ps = short_paths) {
  path = shortest.path.recursive(arm(at_tab[parentid, 10:11]),
                                 arm(at_tab[deadid, 10:11]))
  pathids = rep(0, nrow(path))
  for (i in 1:length(pathids)) {
    pathids[i] = which(at_tab[, 10] == path[i, 1] & at_tab[, 11] == path[i, 2])
  }
  if (length(pathids) > 2) {
    at_tab[pathids[3:length(pathids)], c(2, 3, 5, 7, 8)] =
      at_tab[pathids[2:(length(pathids) - 1)], c(2, 3, 5, 7, 8)]
  }
  daughterid = pathids[2]

  at_tab[daughterid, c(2, 3, 4, 7, 8)] = c(
    draw.hunger.threshold.stdnorm(at_tab[parentid, 3]),
    at_tab[parentid, c(5, 7, 8)])

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

# Hex Plot
plot.hex.net2 = function(double = TRUE, remake = TRUE, at_tab = atr) {

  # code for labeling by id: 'text(atr[, 10], atr[, 11], atr[, 1], cex = .2)'

  # This weird scaling reflects the actual effect of V on K+ uptake
  voltage_colors = floor(abs(at_tab[, 8])*10 - 1519)
  voltage_colors[which(voltage_colors > 81)] = 81
  voltage_colors[which(voltage_colors < 1)] = 1
  frame_color1 = (volt_grad[c(21,101)])[voltage_colors]
  point_color1 = (volt_grad[c(1:81)])[voltage_colors]

  # potassium_colors = floor(at_tab[, 5]) - 190
  # potassium_colors[which(potassium_colors > 140)] = 140
  # potassium_colors[which(potassium_colors < 1)] = 1
  # frame_color2 = (pot_grad[c(21:160)])[potassium_colors]
  # point_color2 = (pot_grad[c(1:140)])[potassium_colors]

  # sig_colors = c("purple", "grey")[(at_tab[, 7]  <= at_tab[, 2]) + 1]
  sig_colors = c("purple", "grey")[is.element(at_tab[, 1], find.sigs.fv()) + 1]

  glutamate_colors = floor(at_tab[, 7]*100) - 20
  glutamate_colors[which(glutamate_colors > 450)] = 450
  glutamate_colors[which(glutamate_colors < 1)] = 1
  frame_color3 = (glu_grad[c(41:490)])[glutamate_colors]
  point_color3 = (glu_grad[c(1:450)])[glutamate_colors]

  threshold_colors = floor(at_tab[, 2]*100) + 30
  threshold_colors[which(threshold_colors <= 30)] = 1
  frame_color4 = (thr_grad[c(31:360)])[threshold_colors]
  point_color4 = (thr_grad[c(1:331)])[threshold_colors]

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
  if (double == FALSE) {
    plot(1, xlim = c(-r, r), ylim = c(-r, r), xlab = "", ylab = "",
         axes = FALSE, frame.plot = FALSE, asp = 1)
    for (i in 1:(length(hex_vectors)/12)) {
      polygon(x = hex_vectors[i, 1:6], y = hex_vectors[i, 7:12],
              border = NA, col = point_color1[i], lty = par("lty"),
              fillOddEven = TRUE)
    }
  } else { # plot signaling and glu
    par(mfrow = c(2, 3), oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 1))

    # Signaling
    plot(1, xlim = c(-r, r), ylim = c(-r, r), xlab = "", ylab = "",
         axes = FALSE, frame.plot = FALSE, asp = 1, main = "Signaling")
    for (i in 1:(length(hex_vectors)/12)) {
      polygon(x = hex_vectors[i, 1:6], y = hex_vectors[i, 7:12],
              border = NA, col = point_color1[i], lty = par("lty"),
              fillOddEven = TRUE)
    }

    # # Potassium
    # plot(1, xlim = c(1, N), ylim = c(M, 0.5), xlab = "", ylab = "",
    #      axes = FALSE, frame.plot = FALSE, asp = 1)
    # for (i in 1:(length(hex_vectors)/12)) {
    #   polygon(x = hex_vectors[i, 1:6], y = hex_vectors[i, 7:12],
    #         border = point_color2[i], col = point_color2[i], lty = par("lty"),
    #         fillOddEven = TRUE)
    # }

    # Signaling
    plot(1, xlim = c(-r, r), ylim = c(-r, r), xlab = "", ylab = "",
         axes = FALSE, frame.plot = FALSE, asp = 1, main = "F(V) Signaling")
    for (i in 1:(length(hex_vectors)/12)) {
      polygon(x = hex_vectors[i, 1:6], y = hex_vectors[i, 7:12],
              border = NA, col = sig_colors[i], lty = par("lty"),
              fillOddEven = TRUE)
    }

    # Glutamate
    plot(1, xlim = c(-r, r), ylim = c(-r, r), xlab = "", ylab = "",
         axes = FALSE, frame.plot = FALSE, asp = 1, main = "Glutamate")
    for (i in 1:(length(hex_vectors)/12)) {
      polygon(x = hex_vectors[i, 1:6], y = hex_vectors[i, 7:12],
              border = point_color3[i], col = point_color3[i], lty = par("lty"),
              fillOddEven = TRUE)
    }

    # # Threshold
    # plot(1, xlim = c(1, N), ylim = c(M, 0.5), xlab = "", ylab = "",
    #      axes = FALSE, frame.plot = FALSE, asp = 1)
    # for (i in 1:(length(hex_vectors)/12)) {
    #   polygon(x = hex_vectors[i, 1:6], y = hex_vectors[i, 7:12],
    #           border = point_color4[i], col = point_color4[i],
    #           lty = par("lty"), fillOddEven = TRUE)
    # }
    # dens = density(Eff(at_tab[, 8]))
    # plot(density(Eff(at_tab[, 8]), from = 0, to = 1), xlim = c(0, 1),
    #      ylab = "Density", xlab = "F(V)")
    # abline(v = dens$x[
    #   which(
    #     is.element(dens$y,
    #                sort(dens$y[which(diff(sign(diff(dens$y))) < 0) + 1],
    #                     decreasing = TRUE)[1:2]
    # ))])

    # Plot F(V) time series for nonsignalers
    plot(x = trackers[1:(stepwise_timer - 1), 18],
         y = trackers[1:(stepwise_timer - 1), 6],
         xlab = "Time", ylab = "Nonsignaler F(V)",
         type = "l", ylim = c(0, 1))
    lines(trackers[1:(stepwise_timer - 1), 3], lty = "dotted")
    # lines(sig_perc_hist, lty = "dashed")
    lines(trackers[1:(stepwise_timer - 1), 7]/3, lty = "dotdash")
    legend("topright", c("F(V)", "Signaling Rate", "% Refractory",
                         "Signaler Threshold (/3)"),
           lty = c("solid", "dotted", "dashed", "dotdash"))

    # Plot F(V) division
    dens = density(Eff(at_tab[, 8]))
    plot(dens, xlab = "F(V)", main = "", xlim = c(0, 1))
    local_maxima <- which(diff(sign(diff(dens$y))) < 0) + 1
    abline(v = dens$x[
      which(
        is.element(dens$y,
                   sort(dens$y[which(diff(sign(diff(dens$y))) < 0) + 1],
                        decreasing = TRUE)[1:2]
        )
      )
    ])

    # plot(density(rowSums(sig_hist)/sig_history_length),
    #      xlab = "Percent Time Signaling")

    # Plot threshold distribution
    plot(density(at_tab[, 2]), xlab = "Threshold", main = "")

    # short_track = trackers %>% as.data.frame() %>% subset(Time != 0)
    # if (length(short_track != 0)) {
    #   # Voltage w time
    #   plot(x = short_track$Time,
    #        y = short_track$Median_V,
    #        type = "l", ylim = c(-156, -150),
    #        ylab = "Membrane Potential", xlab = "Time")
    # }
    par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(5.1, 4.1, 4.1, 2.1))
  }
}

plot.hex.net.quickndirty = function(double = FALSE, remake = FALSE,
                                    as.hex = TRUE, at_tab = atr,
                                    plot.boundary = FALSE) {

  # code for labeling by id: 'text(atr[, 10], atr[, 11], atr[, 1], cex = .2)'

  # This weird scaling reflects the actual effect of V on K+ uptake
  voltage_colors = floor(abs(at_tab[, 8])*10 - 1519)
  voltage_colors[which(voltage_colors > 81)] = 81
  voltage_colors[which(voltage_colors < 1)] = 1
  frame_color1 = (volt_grad[c(21,101)])[voltage_colors]
  point_color1 = (volt_grad[c(1:81)])[voltage_colors]

  ### Remove this lol
  # glu_radius = find.glu.radius()
  # overlay = which(atr[, 9] == glu_radius)
  # frame_color1[overlay] = "yellow"
  # point_color1[overlay] = "yellow"

  if (double) {
    sig_colors = c("purple", "grey")[is.element(at_tab[, 1], find.sigs.fv()) + 1]

    glutamate_colors = floor(at_tab[, 7]*100) - 20
    glutamate_colors[which(glutamate_colors > 450)] = 450
    glutamate_colors[which(glutamate_colors < 1)] = 1
    frame_color3 = (glu_grad[c(41:490)])[glutamate_colors]
    point_color3 = (glu_grad[c(1:450)])[glutamate_colors]

    threshold_colors = floor(at_tab[, 2]*100) + 30
    threshold_colors[which(threshold_colors <= 30)] = 1
    frame_color4 = (thr_grad[c(31:360)])[threshold_colors]
    point_color4 = (thr_grad[c(1:331)])[threshold_colors]
  }

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

  if (double) {
    par(mfrow = c(2, 3), oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 1))

    # Signaling
    plot(1, col = "white", xlim = c(-r, r), ylim = c(-r, r),
         xlab = "", ylab = "", axes = FALSE, frame.plot = FALSE, asp = 1,
         main = "Signaling")
    most_common = most.freq.val(point_color1)
    not_common = which(point_color1 != most_common)
    polygon(x = c(0, -r - 0.64, -r - 0.64, 0, r - 0.3557, r - 0.3557),
            y = c(min(at_tab[, 11]) - 0.8,
                  floor(-r/2) + 0.5,
                  floor(r/2) + 0.5,
                  max(at_tab[, 11]) + 0.8,
                  floor(r/2) + 0.5,
                  floor(-r/2) + 0.5),
            col = most_common, border = NA)
    for (i in not_common) {
      polygon(x = hex_vectors[i, 1:6], y = hex_vectors[i, 7:12],
              border = NA, col = point_color1[i], lty = par("lty"),
              fillOddEven = TRUE)
    }

    # F_v Signaling
    plot(1, col = "white", xlim = c(-r, r), ylim = c(-r, r), xlab = "",
         ylab = "", axes = FALSE, frame.plot = FALSE, asp = 1,
         main = "Signaling")
    most_common = most.freq.val(sig_colors)
    not_common = which(sig_colors != most_common)
    polygon(x = c(0, -r - 0.64, -r - 0.64, 0, r - 0.3557, r - 0.3557),
            y = c(min(at_tab[, 11]) - 0.8,
                  floor(-r/2) + 0.5,
                  floor(r/2) + 0.5,
                  max(at_tab[, 11]) + 0.8,
                  floor(r/2) + 0.5,
                  floor(-r/2) + 0.5),
            col = most_common, border = NA)
    for (i in not_common) {
      polygon(x = hex_vectors[i, 1:6], y = hex_vectors[i, 7:12],
              border = NA, col = sig_colors[i], lty = par("lty"),
              fillOddEven = TRUE)
    }

    # Glutamate
    plot(1, col = "white", xlim = c(-r, r), ylim = c(-r, r), xlab = "",
         ylab = "", axes = FALSE, frame.plot = FALSE, asp = 1,
         main = "Signaling")
    most_common = most.freq.val(point_color3)
    not_common = which(point_color3 != most_common)
    polygon(x = c(0, -r - 0.64, -r - 0.64, 0, r - 0.3557, r - 0.3557),
            y = c(min(at_tab[, 11]) - 0.8,
                  floor(-r/2) + 0.5,
                  floor(r/2) + 0.5,
                  max(at_tab[, 11]) + 0.8,
                  floor(r/2) + 0.5,
                  floor(-r/2) + 0.5),
            col = most_common, border = NA)
    for (i in not_common) {
      polygon(x = hex_vectors[i, 1:6], y = hex_vectors[i, 7:12],
              border = NA, col = point_color3[i], lty = par("lty"),
              fillOddEven = TRUE)
    }

    # Plot F(V) time series for nonsignalers
    plot(x = trackers[1:(stepwise_timer - 1), 18],
         y = trackers[1:(stepwise_timer - 1), 6],
         xlab = "Time", ylab = "Nonsignaler F(V)",
         type = "l", ylim = c(0, 1))
    lines(trackers[1:(stepwise_timer - 1), 3], lty = "dotted")
    # lines(sig_perc_hist, lty = "dashed")
    lines(trackers[1:(stepwise_timer - 1), 7]/3, lty = "dotdash")
    legend("topright", c("F(V)", "Signaling Rate", "% Refractory",
                         "Signaler Threshold (/3)"),
           lty = c("solid", "dotted", "dashed", "dotdash"))

    # Plot F(V) division
    dens = density(Eff(at_tab[, 8]))
    plot(dens, xlab = "F(V)", main = "", xlim = c(0, 1))
    local_maxima <- which(diff(sign(diff(dens$y))) < 0) + 1
    abline(v = dens$x[
      which(
        is.element(dens$y,
                   sort(dens$y[which(diff(sign(diff(dens$y))) < 0) + 1],
                        decreasing = TRUE)[1:2]
        )
      )
    ])

    # Plot threshold distribution
    plot(density(at_tab[, 2]), xlab = "Threshold", main = "")

    par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(5.1, 4.1, 4.1, 2.1))
  } else {
    # Signaling
    point_color1[unoccupied] = "white"
    # most_common = most.freq.val(point_color1)
    # not_common = which(point_color1 != most_common)

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
}

flexy.plot = function(cat, altcols = NA, are.unoccupied = NA, overlay = NaN,
                      overlay_col = "black", at_tab = atr, remake = TRUE,
                      do.frame = FALSE) {
  if (length(altcols) > 1) {
    colors = altcols[cat + 1] # max(cat) + 1]
  } else {
    colors = rainbow(max(cat)*2 + 5)[cat + max(cat) + 1]
  }

  if (length(are.unoccupied) > 1) {
    colors[are.unoccupied] = "white"
  }

  frame_color1 = colors
  frame_color1[which(frame_color1 != "white")] = "black"
  point_color1 = colors

  if (sum(is.nan(overlay)) == 0) {
    frame_color1[overlay] = overlay_col
    point_color1[overlay] = overlay_col
  }

  if (remake == TRUE) {
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
  # most_common = most.freq.val(point_color1)
  # not_common = which(point_color1 != most_common)

  plot(1, col = "white", xlim = c(-r, r), ylim = c(r, -r),
       xlab = "", ylab = "", axes = FALSE, frame.plot = FALSE, asp = 1)

  # polygon(x = c(0, -r - 0.64, -r - 0.64, 0, r + 0.64, r + 0.64),
  #         y = c(min(at_tab[, 11]) - 0.8,
  #               -r/2 - 0.5,
  #               r/2 + 0.5,
  #               max(at_tab[, 11]) + 0.8,
  #               r/2 + 0.5,
  #               -r/2 - 0.5),
  #         col = most_common, border = NA)

  if (do.frame) {
    for (i in 1:length(at_tab[, 1])){ #not_common) {
      polygon(x = hex_vectors[i, 1:6], y = hex_vectors[i, 7:12],
              border = frame_color1[i], lwd = 0.2, col = point_color1[i],
              lty = par("lty"),
              fillOddEven = TRUE)
    }
  } else {
    if (length(are.unoccupied) == 1) {
      for (i in 1:length(at_tab[, 1])){ #not_common) {
        polygon(x = hex_vectors[i, 1:6], y = hex_vectors[i, 7:12],
                border = NA, col = point_color1[i], lty = par("lty"),
                fillOddEven = TRUE)
      }
    } else {
      for (i in which(1:length(at_tab[, 1]) %in% are.unoccupied == FALSE)) {
        polygon(x = hex_vectors[i, 1:6], y = hex_vectors[i, 7:12],
                border = NA, col = point_color1[i], lty = par("lty"),
                fillOddEven = TRUE)
      }
    }
  }
}

# text(at_tab[, 10], at_tab[, 11], at_tab[, 1], cex = .2, col = "yellow")

# This is for plotting as a ggplot for making an animation
# Also a lot faster than my other code and produces a prettier plot
plot.hex.net.gg = function(remake = FALSE, as.hex = TRUE,
                           return.color.df = FALSE, at_tab = atr,
                           cat = NA, altcols = NA, are.unoccupied = NA,
                           plot.boundary = FALSE){

  # code for labeling by id: 'text(atr[, 10], atr[, 11], atr[, 1], cex = .2)'

  # This weird scaling reflects the actual effect of V on K+ uptake
  voltage_colors = floor(abs(at_tab[, 8])*10 - 1519)
  voltage_colors[which(voltage_colors > 81)] = 81
  voltage_colors[which(voltage_colors < 1)] = 1
  frame_color1 = (volt_grad[c(21,101)])[voltage_colors]
  point_color1 = (volt_grad[c(1:81)])[voltage_colors]

  if (length(altcols) > 1) {
      colors = altcols[cat + max(cat) + 1]

    if (length(are.unoccupied) > 1) {
      colors[are.unoccupied] = "#FFFFFF"
      occupied = which(1:length(at_tab[, 1]) %in% are.unoccupied == FALSE)
    }

    frame_color1 = colors
    point_color1 = colors
  }

  if (plot.boundary == TRUE) {
    bounds = which(atr[, 9] %in% (find.glu.radius() + c(-1, 0, 1)))
    frame_color1[bounds] = "yellow"
    point_color1[bounds] = "yellow"
  }

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
    hexdat_x <- reshape2::melt(cbind(id = 1:length(at_tab[, 1]), as.data.frame(hex_x)),
                     id.vars = "id", value.name = "x")
    hexdat_y <- reshape2::melt(cbind(id = 1:length(at_tab[, 1]), as.data.frame(hex_y)),
                     id.vars = "id", value.name = "y")

    # Merge them into the same dataframe
    hexdat <- merge(hexdat_x, hexdat_y)
  }

  colors = cbind(id = at_tab[, 1], col = point_color1)
  plot_df = merge(hexdat, colors)
  plot_df = plot_df[plot_df$id %in% occupied, ]

  if (return.color.df) {
    return(plot_df)
  }

  p1 = ggplot(plot_df, aes(x, y)) +
    geom_polygon(aes(group = id), fill = plot_df$col, colour = "#00000000") +
    coord_fixed(ratio = max(at_tab[occupied, 10])/max(at_tab[occupied, 11])) +
    theme_void()

  return(p1)
}

plot.a.tiff.please = function(title,
                          location = "Plot_Code/Final_Plots/Images/",
                          plot.boundary = FALSE) {
  tiff(paste0(location, title, ".tiff"), width = 5.2, height = 5.2,
       units = "in", res = 600)
  plot.hex.net.quickndirty(plot.boundary = plot.boundary)
  dev.off()
}

plot.singlecell.trajs = function(n = 8) {
  if (stepwise_timer == 1) {
    cell_ids_sct <<- c(sample(find.external(), n - 2),
                       sample(find.internal(), 2))
    k_i_mat <<- matrix(atr[cell_ids_sct, 5], nrow = n, ncol = 1)
    g_i_mat <<- matrix(atr[cell_ids_sct, 7], nrow = n, ncol = 1)
    v_mat <<- matrix(atr[cell_ids_sct, 8], nrow = n, ncol = 1)
  } else {
    k_i_mat <<- cbind(k_i_mat, atr[cell_ids_sct, 5])
    g_i_mat <<- cbind(g_i_mat, atr[cell_ids_sct, 7])
    v_mat <<- cbind(v_mat, atr[cell_ids_sct, 8])
  }

  par(mfrow = c(4, n), mar = c(0.5, 2, 2, 1))

  pre = min(stepwise_timer - 1, 100)

  for (i in 1:n) {
    plot(g_i_mat[i, (stepwise_timer - pre):stepwise_timer],
         type = "l", main = paste("G_i", i),
         xlab = "", xaxt = "n")
  }
  for (i in 1:n) {
    plot(k_i_mat[i, (stepwise_timer - pre):stepwise_timer],
         type = "l", main = paste("K_i", i),
         xlab = "", xaxt = "n")
  }
  for (i in 1:n) {
    plot(v_mat[i, (stepwise_timer - pre):stepwise_timer],
         type = "l", main = paste("V", i),
         xlab = "", xaxt = "n")
  }
  plot(trackers$Out_fv_sig[1:stepwise_timer], ylim = c(0, 0.43),
       type = "l", main = "Out_FV")
  plot(trackers$Out_sig[1:stepwise_timer], type = "l", main = "Out_SIG")
  plot(trackers$In_fv_sig[1:stepwise_timer], type = "l", main = "In_FV")
  plot(trackers$In_sig[1:stepwise_timer], type = "l", main = "In_SIG")
  if (stepwise_timer > 100) {
    plot(trackers$Out_sig_inher[1:stepwise_timer],
         type = "l", main = "Sig_Inher")
    plot(trackers$Out_non_inher[1:stepwise_timer],
         type = "l", main = "Non_Inher")
  }

  par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
}

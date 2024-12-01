grow.a.biofilm3 = function(method = c("entire", "first_layer", "add_layer")) {
  if (method == "entire" | method == "first_layer") {
    atr = make.hex.atr()
    atr[2:length(atr[, 1]), c(2,3)] = NA
    atr[1, 2:3] = mean(c(upper_thresh, lower_thresh))
    occupied = 1

    empty_neigh_vec = apply(
      atr, 1,
      function(row) {sum(is.na(row[12:17]) == FALSE)}
    )

    empty_neigh_vec[2:7] %-=% 1

    parent_matrix = matrix(c(1, 1), nrow = 1, ncol = 2)
    colnames(parent_matrix) = c("Offspring", "Parent")
  }
  if (method == "entire") {
    while (length(occupied)/length(atr[, 1]) < prop_occupied) {
      border_cells = get.border.cells()

      parent = sample(border_cells, 1,
                      prob = 1/(atr[border_cells, 9] + 1))

      potential_daughter = atr[parent, 12:17][
        which(is.element(atr[parent, 12:17], occupied) == FALSE &
                is.na(atr[parent, 12:17]) == FALSE)]
      daughter = potential_daughter[sample(length(potential_daughter), 1)]

      parent_matrix = rbind(parent_matrix, c(daughter, parent))

      atr[daughter, c(2, 3)] = draw.hunger.threshold.stdnorm(atr[parent, 3])

      occupied = c(occupied, daughter)

      empty_neigh_vec[atr[daughter, 12:17][which(is.element(
        atr[daughter, 12:17], occupied))]] %-=% 1

      if (length(unique(occupied)) != length(occupied)) {break()}

      if (length(occupied) %% 1000 == 1) {
        # flexy.plot(is.element(atr[, 1], occupied), at_tab = atr)
        print(length(occupied))
      }
    }
  } else if (method == "add_layer") {
    border_cells = get.border.cells()

    new_baby_num  = ceiling(length(border_cells)*gens_per_tick)

    # identify parents
    parents = sample(border_cells, new_baby_num,
                     prob = 1/(atr[border_cells, 9] + 1))

    # draw daughters efficiently with apply
    if (length(parents) == 1) {
      daughters = sample.elig.daughter(atr[parents, ])
    } else {
      daughters = apply(atr[parents, ], 1, sample.elig.daughter)
      while (length(daughters) != length(unique(daughters))) {
        parents = sample(border_cells, new_baby_num,
                         prob = 1/(atr[border_cells, 9] + 1))
        daughters = apply(atr[parents, ], 1, sample.elig.daughter)
      }
    }

    print(daughters)
    # remove duplicate daughters

    parent_matrix = rbind(parent_matrix, cbind(daughters, parents))

    atr[daughters, c(2, 3, 5, 7, 8)] = cbind(
      draw.hunger.threshold.stdnorm(atr[parents, 3]),
      rep(300, new_baby_num),
      sample(seq(upper_thresh, G_m, by = 0.01), new_baby_num),
      rep(V_0, new_baby_num)
    )

    occupied = c(occupied, daughters)

    daughter_neighs = c(atr[daughters, 12:17])

    # account for new daughters in empty_neigh_vec
    dnc = table(daughter_neighs)

    empty_neigh_vec[unique(daughter_neighs)] =
      empty_neigh_vec[unique(daughter_neighs)] - dnc

    if (length(unique(occupied)) != length(occupied)) {
      warning("Occupied duplicates")
      print(occupied)
      break()
    }

  }
  atr[which(is.element(atr[, 1], occupied) == FALSE), c(5, 7, 8)] = 0
  atr[which(is.element(atr[, 1], occupied) == FALSE), 2:3] = -1

  parent_matrix <<- parent_matrix
  empty_neigh_vec <<- empty_neigh_vec
  occupied <<- occupied
  return(atr)
}


sample.elig.daughter = function(row) {
  index_list = which(is.element(row[12:17], occupied) == FALSE &
                       is.na(row[12:17]) == FALSE)
  if (length(index_list) == 1) {
    return(row[11 + index_list])
  } else {
    return(row[11 + sample(index_list, 1)])
  }
}

count.empty.neighs = function(row) {
  sum(is.element(row[12:17], occupied) == FALSE & is.na(row[12:17]) == FALSE)
}

get.border.cells = function() {
  intersect(which(empty_neigh_vec > 0), occupied)
}


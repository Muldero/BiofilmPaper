setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path), ".."))

temp.sig.def = function(vec, subset = 1:length(vec)) {
  F_v = Eff(vec[occupied])
  mid_val = weighted.mean(c(mean(F_v), 1), c(2/3, 1/3))

  # Now identify the maxima for cells below the midpoint
  lower_half_fv = F_v[which(F_v < mid_val)]
  dens = density(lower_half_fv, adjust = 0.2)
  xval = dens$x[which.max(dens$y)]

  # Now define our cutoff to be a weighted average between x and midval
  cutoff = weighted.mean(c(xval, 1), c(3/4, 1/4))

  return(subset %in% intersect(occupied[which(F_v > cutoff)], subset))
}

return_folder = "Output"
tracker_replications = c(1:20)

return_table = matrix(ncol = 6, nrow = 11)

colnames(return_table) = c("Outer", "Inner", "Total",
                           "O_error", "In_error", "Tot_error")
rownames(return_table) = c(
  "Signaling Fraction", "Pairwise Signaler Recurrence",
  "Pairwise Non-signaler Recurrence", "Total Signaler Recurrence",
  "Total Non-signaler Recurrence", "Pairwise Consistent Signaling Fraction",
  "Pairwise Consistent Non-signaling Fraction",
  "Pairwise Inconsistent Fraction", "Total Consistent Signaling Fraction",
  "Total Consistent Non-signaling Fraction", "Total Inconsistent Fraction"
)

################################################################################
### Data from Trackers #########################################################
################################################################################

trackers = read.csv(paste0("Output/POST_TRACK/Post_track_", 1,
                           ".csv"), header = TRUE)

store_o = matrix(ncol = 0, nrow = 3)
store_i = matrix(ncol = 0, nrow = 3)
store_t = matrix(ncol = 0, nrow = 3)
rownames(store_o) = colnames(trackers)[c(16, 26, 27)]
rownames(store_i) = colnames(trackers)[c(4, 14, 15)]
rownames(store_t) = colnames(trackers)[c(28, 38, 39)]

for (i in tracker_replications) {
  trackers = read.csv(paste0("Output/POST_TRACK/Post_track_", i,
                             ".csv"), header = TRUE)

  po = get.only.max.peaks(trackers$Out_fv_sig[1501:3000]) + 1500
  pi = get.only.max.peaks(trackers$In_fv_sig[1501:3000],
                          midline = 0.47) + 1500
  pt = get.only.max.peaks(trackers$Tot_fv_sig[1501:3000],
                          midline = 0.2) + 1500

  store_o = cbind(store_o, t(trackers[po, c(16, 26, 27)]))
  store_i = cbind(store_i, t(trackers[pi, c(4, 14, 15)]))
  store_t = cbind(store_t, t(trackers[pt, c(28, 38, 39)]))
}

return_table[1:3, c(1, 4)] = c(
  rowMeans(store_o),
  sd(store_o[1, ]),
  sd(store_o[2, ]),
  sd(store_o[3, ])
)

return_table[1:3, c(2, 5)] = c(
  rowMeans(store_i),
  sd(store_i[1, ]),
  sd(store_i[2, ]),
  sd(store_i[3, ])
)

return_table[1:3, c(3, 6)] = c(
  rowMeans(store_t),
  sd(store_t[1, ]),
  sd(store_t[2, ]),
  sd(store_t[3, ])
)

################################################################################
### Data from Membrane Potential Matrix ########################################
################################################################################

outer_storage = matrix(ncol = 5, nrow = 27)
rownames(outer_storage) = c(
  "Counts_Out", "Counts_In", "Counts_Tot", # 1:3

  "Pair_cs_out", "Pair_cn_out", "Pair_cns_out", # 4:6
  "Pair_cs_in", "Pair_cn_in", "Pair_cns_in", # 7:9
  "Pair_cs_tot", "Pair_cn_tot", "Pair_cns_tot", # 10:12

  "Tot_cs_out", "Tot_cn_out", "Tot_cns_out", # 13:15
  "Tot_cs_in", "Tot_cn_in", "Tot_cns_in", # 16:18
  "Tot_cs_tot", "Tot_cn_tot", "Tot_cns_tot", # 19:21

  "Tot_rs_out", "Tot_rn_out", # 22:23
  "Tot_rs_in", "Tot_rn_in", # 24:25
  "Tot_rs_tot", "Tot_rn_tot" # 26:27
)

for (i in 1:5) {
  sig_traj = readRDS(paste0(
    "Output/Other/Membrane_Potential_Matrix_",
    i, ".RDS"))[1501:3000, ]

  external = readRDS(paste0(
    "Output/Other/Membrane_Potential_External_", i, ".RDS"))
  sig_traj_ex = sig_traj[, external]

  occupied = which(colMeans(sig_traj) != 0)
  internal = which(!(1:dim(sig_traj)[2] %in% external))
  internal = intersect(internal, occupied)

  trackers = readRDS(paste0(
    "Output/Other/Membrane_Potential_Tracker_", i, ".RDS"))
  trackers = filter(trackers, trackers$Time > 1500)

  parent_matrix = readRDS(paste0(
    "Output/Other/Membrane_Potential_Parent_Matrix_",
    i, ".RDS"))

  peaks = get.only.max.peaks(trackers$Out_fv_sig)
  dips = get.only.max.peaks(trackers$Out_fv_sig, midline = 0.25, top = FALSE)
  num_osc = length(dips)

  peaks_in = get.only.max.peaks(trackers$In_fv_sig, midline = 0.47)
  dips_in = get.only.max.peaks(trackers$In_fv_sig, midline = 0.45, top = FALSE)
  num_osc_in = length(dips_in)

  peaks_tot = get.only.max.peaks(trackers$Tot_fv_sig, midline = 0.2)
  dips_tot = get.only.max.peaks(trackers$Tot_fv_sig, midline = 0.25, top = FALSE)
  num_osc_tot = length(dips_tot)

  ##############################################################################
  ### Pairwise Consistency
  ##############################################################################

  ### Out

  ss_count = 0
  nn_count = 0
  ns_count = 0

  sigs = temp.sig.def(sig_traj[peaks[1], ], subset = external)

  for (j in 2:length(peaks)) {
    prev_sigs = sigs
    sigs = temp.sig.def(sig_traj[peaks[j], ], subset = external)

    ss_count = ss_count + sum((sigs + prev_sigs) == 2)
    ns_count = ns_count + sum((sigs + prev_sigs) == 1)
    nn_count = nn_count + sum((sigs + prev_sigs) == 0)
  }

  outer_storage[1, i] = length(external)
  outer_storage[4:6, i] = c(ss_count, nn_count, ns_count)

  ### In

  ss_count = 0
  nn_count = 0
  ns_count = 0

  sigs = temp.sig.def(sig_traj[peaks_in[1], ], subset = internal)

  for (j in 2:length(peaks_in)) {
    prev_sigs = sigs
    sigs = temp.sig.def(sig_traj[peaks_in[j], ], subset = internal)

    ss_count = ss_count + sum((sigs + prev_sigs) == 2)
    ns_count = ns_count + sum((sigs + prev_sigs) == 1)
    nn_count = nn_count + sum((sigs + prev_sigs) == 0)
  }

  outer_storage[2, i] = length(internal)
  outer_storage[7:9, i] = c(ss_count, nn_count, ns_count)

  ### Tot

  ss_count = 0
  nn_count = 0
  ns_count = 0

  sigs = temp.sig.def(sig_traj[peaks_tot[1], ], subset = occupied)

  for (j in 2:length(peaks_tot)) {
    prev_sigs = sigs
    sigs = temp.sig.def(sig_traj[peaks_tot[j], ], subset = occupied)

    ss_count = ss_count + sum((sigs + prev_sigs) == 2)
    ns_count = ns_count + sum((sigs + prev_sigs) == 1)
    nn_count = nn_count + sum((sigs + prev_sigs) == 0)
  }

  outer_storage[3, i] = length(occupied)
  outer_storage[10:12, i] = c(ss_count, nn_count, ns_count)


  ##############################################################################
  ### Total Consistency
  ##############################################################################

  ### Out

  Long_Sig_Out = 1:length(external)

  for (j in 1:length(external)) {
    k = external[j]
    traj_peaks = which(sig_traj[, k] < -390)

    osc_count = 0

    for (x in 1:(num_osc - 1)) {
      if (sum(traj_peaks > dips[x] & traj_peaks < dips[x + 1]) > 0) {
        osc_count = osc_count + 1
      }
    }

    Long_Sig_Out[j] = osc_count/(num_osc - 1)

    if(Long_Sig_Out[j] > 1) {
      Long_Sig_Out[j] = 1
    }
  }

  outer_storage[13:15, i] = c(sum(Long_Sig_Out >= 0.9),
                              sum(Long_Sig_Out <= 0.1),
                              sum(Long_Sig_Out < 0.9 & Long_Sig_Out > 0.1))

  ### In

  Long_Sig_In = 1:length(internal)

  for (j in 1:length(internal)) {
    k = internal[j]
    traj_peaks = which(sig_traj[, k] < -230)

    osc_count = 0

    for (x in 1:(num_osc_in - 1)) {
      if (sum(traj_peaks > dips_in[x] & traj_peaks < dips_in[x + 1]) > 0) {
        osc_count = osc_count + 1
      }
    }

    Long_Sig_In[j] = osc_count/(num_osc_in - 1)

    if(Long_Sig_In[j] > 1) {
      Long_Sig_In[j] = 1
    }
  }

  outer_storage[16:18, i] = c(sum(Long_Sig_In >= 0.9), sum(Long_Sig_In <= 0.1),
                              sum(Long_Sig_In < 0.9 & Long_Sig_In > 0.1))

  ### Tot

  Long_Sig_Tot = 1:length(occupied)

  for (j in 1:length(occupied)) {
    k = occupied[j]
    traj_peaks = which(sig_traj[, k] < -230)

    osc_count = 0

    for (x in 1:(num_osc_tot - 1)) {
      if (sum(traj_peaks > dips_tot[x] & traj_peaks < dips_tot[x + 1]) > 0) {
        osc_count = osc_count + 1
      }
    }

    Long_Sig_Tot[j] = osc_count/(num_osc_tot - 1)

    if(Long_Sig_Tot[j] > 1) {
      Long_Sig_Tot[j] = 1
    }
  }

  outer_storage[19:21, i] = c(sum(Long_Sig_Tot >= 0.9),
                              sum(Long_Sig_Tot <= 0.1),
                              sum(Long_Sig_Tot < 0.9 & Long_Sig_Tot > 0.1))

  ##############################################################################
  ### Total Recurrence
  ##############################################################################

  ### Out

  ids = external
  out_sigs = which(Long_Sig_Out > 0.9)
  out_non = which(Long_Sig_Out < 0.1)
  ids = ids[union(out_sigs, out_non)]

  parent_matrix_t = parent_matrix[which(is.element(parent_matrix[, 1],
                                                   ids)), ]
  parent_matrix_t = parent_matrix_t[which(is.element(parent_matrix_t[, 2],
                                                     ids)), ]

  daughter_states = is.element(parent_matrix_t[, 1], out_sigs)
  parent_states = is.element(parent_matrix_t[, 2], out_sigs)

  sig_inher_rate = sum(daughter_states == parent_states &
                         daughter_states == 1)/sum(daughter_states == 1)
  non_inher_rate = sum(daughter_states == parent_states &
                         daughter_states == 0)/sum(daughter_states == 0)

  outer_storage[22:23, i] = c(sig_inher_rate, non_inher_rate)

  ### In

  ids = internal
  in_sigs = which(Long_Sig_In > 0.9)
  in_non = which(Long_Sig_In < 0.1)
  ids = ids[union(in_sigs, in_non)]

  parent_matrix_t = parent_matrix[which(is.element(parent_matrix[, 1],
                                                   ids)), ]
  parent_matrix_t = parent_matrix_t[which(is.element(parent_matrix_t[, 2],
                                                     ids)), ]

  daughter_states = is.element(parent_matrix_t[, 1], in_sigs)
  parent_states = is.element(parent_matrix_t[, 2], in_sigs)

  sig_inher_rate = sum(daughter_states == parent_states &
                         daughter_states == 1)/sum(daughter_states == 1)
  non_inher_rate = sum(daughter_states == parent_states &
                         daughter_states == 0)/sum(daughter_states == 0)

  outer_storage[24:25, i] = c(sig_inher_rate, non_inher_rate)

  ### Tot

  ids = occupied
  tot_sigs = which(Long_Sig_Tot > 0.9)
  tot_non = which(Long_Sig_Tot < 0.1)
  ids = ids[union(tot_sigs, tot_non)]

  parent_matrix_t = parent_matrix[which(is.element(parent_matrix[, 1],
                                                   ids)), ]
  parent_matrix_t = parent_matrix_t[which(is.element(parent_matrix_t[, 2],
                                                     ids)), ]

  daughter_states = is.element(parent_matrix_t[, 1], tot_sigs)
  parent_states = is.element(parent_matrix_t[, 2], tot_sigs)

  sig_inher_rate = sum(daughter_states == parent_states &
                         daughter_states == 1)/sum(daughter_states == 1)
  non_inher_rate = sum(daughter_states == parent_states &
                         daughter_states == 0)/sum(daughter_states == 0)

  outer_storage[26:27, i] = c(sig_inher_rate, non_inher_rate)

  print(i)
}

################################################################################
### Storing Everything #########################################################
################################################################################

# Total Signaler Recurrence
return_table[4, 1] = weighted.mean(outer_storage[22, ], outer_storage[1, ])
return_table[4, 2] = weighted.mean(outer_storage[24, ], outer_storage[2, ])
return_table[4, 3] = weighted.mean(outer_storage[26, ], outer_storage[3, ])

# Total Non-signaler Recurrence
return_table[5, 1] = weighted.mean(outer_storage[23, ], outer_storage[1, ])
return_table[5, 2] = weighted.mean(outer_storage[25, ], outer_storage[2, ])
return_table[5, 3] = weighted.mean(outer_storage[27, ], outer_storage[3, ])

# Pairwise Consistent Out
return_table[6:8, 1] = c(sum(outer_storage[4, ]), sum(outer_storage[5, ]),
                         sum(outer_storage[6, ]))/sum(outer_storage[4:6, ])

# Pairwise Consistent In
return_table[6:8, 2] = c(sum(outer_storage[7, ]), sum(outer_storage[8, ]),
                         sum(outer_storage[9, ]))/sum(outer_storage[7:9, ])

# Pairwise Consistent Tot
return_table[6:8, 3] = c(sum(outer_storage[10, ]), sum(outer_storage[11, ]),
                         sum(outer_storage[12, ]))/sum(outer_storage[10:12, ])

# Total Consistent Out
return_table[9:11, 1] = c(sum(outer_storage[13, ]), sum(outer_storage[14, ]),
                         sum(outer_storage[15, ]))/sum(outer_storage[1, ])

# Total Consistent In
return_table[9:11, 2] = c(sum(outer_storage[16, ]), sum(outer_storage[17, ]),
                          sum(outer_storage[18, ]))/sum(outer_storage[2, ])

# Total Consistent Tot
return_table[9:11, 3] = c(sum(outer_storage[19, ]), sum(outer_storage[20, ]),
                          sum(outer_storage[21, ]))/sum(outer_storage[3, ])

# Add in Errors

se = function(p, n) {
  sqrt((
    sum(n*(p - weighted.mean(p, n))^2)/sum(n)
  )/length(p))
}

for (i in 1:3) { # Out, in, total
  # Total Recurrence
  return_table[4, i + 3] = se(outer_storage[20 + 2*i, ], outer_storage[i, ])
  return_table[5, i + 3] = se(outer_storage[21 + 2*i, ], outer_storage[i, ])

  # Pairwise Consistent
  for (j in 6:8) {
    return_table[j, i + 3] = se(outer_storage[(j - 5) + i*3, ] /
                                  sum(outer_storage[i*3 + 1:3, ]),
                                outer_storage[i, ])
  }

  # Total Consistent
  for (j in 9:11) {
    return_table[j, i + 3] = se(outer_storage[(j + 1) + i*3, ] /
                                  sum(outer_storage[i, ]),
                                outer_storage[i, ] * 0.9)
    # multiply by 0.9 to account for not quite all here
  }
}

# Report sd

for (i in 1:3) { # Out, in, total
  # Total Recurrence
  return_table[4, i + 3] = sd(outer_storage[20 + 2*i, ])
  return_table[5, i + 3] = sd(outer_storage[21 + 2*i, ])

  # Pairwise Consistent
  for (j in 6:8) {
    return_table[j, i + 3] = sd(outer_storage[(j - 5) + i*3, ] /
                                  sum(outer_storage[i*3 + 1:3, ]))
  }

  # Total Consistent
  for (j in 9:11) {
    return_table[j, i + 3] = sd(outer_storage[(j + 1) + i*3, ] /
                                  sum(outer_storage[i, ]))
  }
}

return_table = return_table %>% round(3)

View(return_table)

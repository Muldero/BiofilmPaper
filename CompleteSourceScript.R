setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path)))

### Library ####################################################################

source("Packages.R", verbose = FALSE)

options(max.print = 3000, digits = 8)
options(scipen = 8)

plot.gif = FALSE # plot a gif of the model
grow.all.the.way = FALSE # grow the entire biofilm (remove the growht phase)

# memory intensive (!) diagnostic trackers
track.min.glu = FALSE
track.all.sigs = FALSE

ticks = 3000 # how many ticks to run the model for

tpp = 100 # ticks per period (of 2 hours) so divide by tpp/2 for time-dependence

gens_per_tick = 1/40 # number of generations (amount of growth) per tick

source("Functions.R", verbose = FALSE)

return_folder = paste0("Output/")

### Variables ##################################################################

r = 150 # radius of network
init_radius = 50
prop_occupied = 0.75 # proportion of the network that should be occupied if grow

# Observed (true) parameters
K_shock = K_shock_default = -10
K_m = 8 # basal potassium concentration in media
G_m = 30 # basal glutamate concentration in media
V_0 = -150 # base membrane potential

# Potassium parameters
F_param_default = 0.05/3.1 # default 0.05
D_p_default = 0.12*1.5 # 1.5 # default 0.12

G_m_vec = rep(G_m, r*6 + 6)

# Parameters (drawn from Martinez-Corral (2019))
alpha_g = 24 # glu uptake constant
delta_g = 4.8 # glu degradation constant
g_k = 70

# Threshold Parameters
upper_thresh_default = 3
lower_thresh_default = 0
thresh_herit = 1
sd2_default = 1 # standard deviation of genotype

# Range of threshold parameters for alternate runs
ut_vec = c(3,    3,   3,   3, 3, 3,   0, 1, 2, 3, 2.1,   1, 2, 2, 3, 3, 5)
lt_vec = c(0,    0,   0,   0, 0, 0,   0, 1, 2, 3, 1.9,   0, 0, 1, 1, 2, 0)
sd_vec = c(0, 0.01, 0.1, 0.5, 1, 2,   0, 0, 0, 0, 0.001,   3, 3, 3, 3, 3, 3)

# range of glutamate parameters for alternate runs
Fp_vec = c(  2,   2,   2,   2,   1.9, 1.7, 2.1, 2.3,   1.8, 1.8, 2.2, 2.2)
Dp_vec = c(1.9, 1.7, 2.1, 2.3,     2,   2,   2,   2,   2.2, 1.8, 2.2, 1.8)

run_vec = c(
  1:20,
  paste0("Thresh[", lt_vec, ",", ut_vec, "]_SD", sd_vec), # 21:37
  paste0("F", Fp_vec, ",", "D_p", Dp_vec), # 38:49
  "KShock_750", # 50
  "LowGlu", # 51
  "Radius250" # 52
)

for (run_count in c(1:52)) {
  # Collect complete data on the following runs:
  if (run_count %in% c(1:5, 27, 29, 32:35)) {
    track.min.glu = TRUE
    track.all.sigs = TRUE
  } else {
    track.min.glu = FALSE
    track.all.sigs = FALSE
  }
  # Normal runs
  if (run_count <= 20) {
    upper_thresh = upper_thresh_default
    lower_thresh = lower_thresh_default
    sd2 = sd2_default
    F_param = F_param_default
    D_p = D_p_default
    K_shock = K_shock_default

    K_m = 8 # basal potassium concentration in media
    G_m = 30 # basal glutamate concentration in media
    G_m_vec = rep(G_m, r*6 + 6)

    source("Master_Model_SourceOnly.R")

  # Alternate threshold runs
  } else if (run_count <= 37) {
    upper_thresh = ut_vec[run_count - 20]
    lower_thresh = lt_vec[run_count - 20]
    sd2 = sd_vec[run_count - 20]
    K_shock = K_shock_default

    F_param = F_param_default
    D_p = D_p_default

    K_m = 8 # basal potassium concentration in media
    G_m = 30 # basal glutamate concentration in media
    G_m_vec = rep(G_m, r*6 + 6)

    source("Master_Model_SourceOnly.R")

  # Alternate Potassium metabolism runs
  } else if (run_count <= 49) {
    upper_thresh = upper_thresh_default
    lower_thresh = lower_thresh_default
    sd2 = sd2_default
    K_shock = K_shock_default

    F_param = 0.05/Fp_vec[run_count - 37]
    D_p = 0.12*Dp_vec[run_count - 37]

    K_m = 8 # basal potassium concentration in media
    G_m = 30 # basal glutamate concentration in media
    G_m_vec = rep(G_m, r*6 + 6)

    source("Master_Model_SourceOnly.R")

  # Depolarize at Tick 750
  } else if (run_vec[run_count] == "KShock_750") {
    upper_thresh = upper_thresh_default
    lower_thresh = lower_thresh_default
    sd2 = sd2_default
    F_param = F_param_default
    D_p = D_p_default

    K_shock = 750
    kshock_level = 300

    K_m = 8 # basal potassium concentration in media
    G_m = 30 # basal glutamate concentration in media
    G_m_vec = rep(G_m, r*6 + 6)

    source("Master_Model_SourceOnly.R")

  # Grow in low glutamate
  } else if (run_vec[run_count] == "LowGlu") {
    upper_thresh = upper_thresh_default
    lower_thresh = lower_thresh_default
    sd2 = sd2_default
    F_param = F_param_default
    D_p = D_p_default

    r = 100 # 150 # radius of network
    init_radius = 25
    prop_occupied = 0.75

    K_shock = K_shock_default
    G_m = 20
    G_m_vec = rep(G_m, r*6 + 6)

    source("Master_Model_SourceOnly.R")

  # Grow to a larger radius
  } else if (run_vec[run_count] == "Radius250") {
    upper_thresh = upper_thresh_default
    lower_thresh = lower_thresh_default
    sd2 = sd2_default
    F_param = F_param_default
    D_p = D_p_default

    r = 250 # radius of network
    init_radius = 75
    prop_occupied = 0.75

    K_shock = K_shock_default
    G_m = 35
    G_m_vec = rep(G_m, r*6 + 6)

    source("Master_Model_SourceOnly.R")
  }

  # Record Output
  write.csv(atr, paste0(return_folder, "POST_ATR/Post_atr_",
                        run_vec[run_count], ".csv"))

  write.csv(who_is_signaling[which(occupied %in% find.external()), ],
            paste0(return_folder, "SIG_TRACKER/SIG_TRACKER_",
                   run_vec[run_count], ".csv"))

  write.csv(trackers, paste0(return_folder, "POST_TRACK/Post_track_",
                             run_vec[run_count], ".csv"))

  if (run_count %in% c(1:5, 27, 29, 32:35)) {
    write.csv(pot_traj, paste0(return_folder, "Other/Potassium_Trajectories",
                               run_vec[run_count], ".csv"), row.names = FALSE)

        write.csv(glu_traj, paste0(return_folder,
                               "Other/Glutamate_Trajectories_",
                               run_vec[run_count], ".csv"), row.names = FALSE)

    write.csv(mem_pot_traj,
              paste0(return_folder, "Other/Membrane_Potential_Trajectories_",
                     run_vec[run_count], ".csv"), row.names = FALSE)

    write.csv(who_is_signaling, paste0(return_folder,
                     "SIG_TRACKER/SIG_TRACKER_TOTAL_", run_vec[run_count],
                     ".csv"), row.names = FALSE)

    write.csv(who_is_signaling[which((occupied %in% find.external()) == 0), ],
              paste0(return_folder, "SIG_TRACKER/SIG_TRACKER_INNER",
                     run_vec[run_count], ".csv"), row.names = FALSE)

    write.csv(min_glu, paste0(return_folder,
                              "Other/Minimum_Glu_Signaling_",
                              run_vec[run_count], ".csv"), row.names = FALSE)

    saveRDS(sig_tracking, paste0(return_folder,
                                 "Other/Membrane_Potential_Matrix_",
                                 run_vec[run_count], ".RDS"))

    saveRDS(trackers, paste0(return_folder,
                             "Other/Membrane_Potential_Tracker_",
                             run_vec[run_count], ".RDS"))

    saveRDS(parent_matrix, paste0(return_folder,
                                  "Other/Membrane_Potential_Parent_Matrix_",
                                  run_vec[run_count], ".RDS"))

    saveRDS(find.external(), paste0(return_folder,
                                    "Other/Membrane_Potential_External_",
                                    run_vec[run_count], ".RDS"))
  } else {
    write.csv(signal_trajectories,
              paste0(return_folder, "SIG_SAMPLE_TRAJ/SIG_SAMPLE_TRAJ_",
                     run_vec[run_count], ".csv"))
  }

  print(run_vec[run_count])
}

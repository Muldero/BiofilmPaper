setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path), ".."))

# This script must be run before plotting code to load necessary packages and
# functions

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

# rows of the atr representing empty hexes
unoccupied = which(atr[, 2] == -1)

# rows of the atr with cells in them
occupied = which(is.element(atr[, 1], unoccupied) == FALSE)

# boundary between interior and exterior
glu_radius = find.glu.radius()
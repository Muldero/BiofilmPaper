list_of_packages <- c(
  "roperators", # allows %+=% operator
  "magrittr", # allows use of pipes, probably not needed
  "extraDistr", # fancy distributions
  "ape", # Newick tree, maybe not needed
  "igraph", # Network plotting
  "plotly", # more plotting stuff
  "latex2exp", # nice graph labels
  "e1071", # Floyd-Warshall algorithm
  "magick", # gif creation
  "ggpubr", # ggplot2 helper
  "boot", # absolutely no clue what I'm using this for
  "progress", # progress bar
  "zoo", # rolling mean function
  "dplyr",
  "tidyverse",
  "reshape2",
  "patchwork",
  "psych", # harmonic mean
  "grid",
  "ggtext",
  "gridExtra", # arrange plots in ggplot
  "data.table"
)

new_packages <- list_of_packages[!(list_of_packages %in%
                                     installed.packages())]

if (length(new_packages) > 0) {install.packages(new_packages)}

lapply(list_of_packages, library, character.only = TRUE)



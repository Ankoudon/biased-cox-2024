################################################################################

# Thesis: "Impact of Correlation of Exposure between Two Consecutive
# Periods On Estimating Vaccine Efficacy", Hiroyasu Ando, Akihiro Nishi
# Purpose: Create network structures

################################################################################

# (1) Clean the environmennt. 
rm(list = ls())

################################################################################

# (2) Read libraries
library(igraph)
library(tidyverse)
library(Matrix)
library(arrow)

################################################################################

# (3) Parameters 

# Number of vertices
people_n <- 10000000
# Number of simulations
session <- 5
# Infectious period (mean) 
i_prob <- 3
# Symptomatic ratio
symptomatic_rate <- 0.55

################################################################################

# (4) Create network structures

# Set the directory
setwd("~/Desktop/biased-cox-2024/script")

for (h in 1:session) {
  
  print(h)
  
  set.seed(h)
  
  # Barabási–Albert model (power = 0.01, m = 2)
  g <- sample_pa(n = people_n, power = 0.65, m = 3, directed = FALSE)  
  # Create adjacency matrix 
  ydata0 <- as_adjacency_matrix(g, type = "both",
                                sparse = igraph_opt("sparsematrices"))
  
  setwd("../mtx") 
  # Save the adjacency matrix
  writeMM(obj = ydata0, file = paste0("h", h, ".mtx"))
  
  ndata <- tibble(
    # Number of vertices
    id = as.integer(1:people_n),
    # People's states, {s, e, i_asym, i_sym_pre, i_sym_post, r_asym, r_sym}
    state = factor("s", levels = c("s",
                                   "e",
                                   "i_asym",
                                   "i_sym_pre",
                                   "i_sym_post",
                                   "r_asym",
                                   "r_sym")),
    # 1-210 round (equivalent to 1-210 days)   
    round = 0,
    # Flag for the new E stage
    new_e = 0,
    # Flag for the new I stage
    new_i = 0,
    # Flag for the new R stage
    new_r = 0,
    # Number of contacts
    contacts = 0,
    # Number of degrees (neighbors)
    degree = degree(g),
    # Session (simulation) number
    session = 0,
    # Vaccine effect, {0: ineffective, 1: effective}
    vaccine_effect = 0,
    # Vaccine status, {-1: not vaccinated, 1: real vaccine, 0: placebo}
    vaccine_status = -1,
    # Date when the vaccine gets effective
    effect_date = 0,
    # Number of contacts when individuals get infected
    required_contacts = 0,
    # Flag for survival analysis, {0: not included, 1: included}
    test_sample = 0,
    # PCR test before the second dose, {0: not pass, 1: pass}
    test_dose2 = 0,
    # Date for the second dose
    dose2_date = 0,
    # Date when participants are put in the analysis sample
    sample_date = 0,
    # Survival length
    survd = 0,
    # State end date
    state_end = 0,
    # True vaccine effect 
    true_effect = 0,
    # Date participants get vaccinated 
    vaccine_date = 0,
    # Assign I period 
    i_period = rgeom(people_n, prob = (1 / i_prob)) + 1,
    # Assign symptomatic and asymptomatic, {0: asymptomatic, 1: symptomatic}
    symptomatic = sample(c(0, 1),
                         people_n,
                         replace = T,
                         prob = c(1 - symptomatic_rate, symptomatic_rate)),
    # Assign pre-symptomatic period
    pre_sym_period = rbinom(people_n,
                            size = i_period,
                            prob = 0.5),
    # Total number of infectious contacts
    inf_exposure = 0)
  
  setwd("../ndata")
  # Save the population data as a parquet file
  write_parquet(as_arrow_table(ndata),
                paste0("h", h, "_ndata.parquet"))
  
  NULL
  
}
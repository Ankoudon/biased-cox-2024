################################################################################

# Thesis: "Impact of Correlation of Exposure between Two Consecutive
# Periods On Estimating Vaccine Efficacy", Hiroyasu Ando, Akihiro Nishi
# Purpose: To simulate the clinical trial data

################################################################################

# (1) Clean the environmenmt. 
rm(list = ls())

################################################################################

# (2) Read libraries

library(igraph)
library(tidyverse)
library(Matrix)
library(arrow)
library(survival)
library(coxphf)

################################################################################

# (3) Parametrers

# Number of simulations
session <- 5
# Number of people
people_n <- 3000000
# Per-contact transmissibility
per_trans <- 0.064
# E period
e_period <- 3
# R period 
r_period <- 0
# Number of rounds (equivalent to days) 
round_num <- 210
# Number of participants 
subjects_n <- 300000
# Length of time for vaccine to be effective  
sample_period <- 42
# Vaccine efficacy 
v_benefit <- 0.941
# External infection per week
external_inf <- 300
# Initial E stage people 
initial_e <- 800
# Initial R stage people 
initial_r = 18262

################################################################################

# (4) Run the clinical trial simulations

# Set the directory
setwd("~/Desktop/biased-cox-2024/script")

################################################################################

for (h in 1:session) {
  
  setwd("~/Desktop/biased-cox-2024/ndata")
  
  set.seed(h)
  
  # Summary of the simulation
  result <- tibble(
    #Number of simulations
    session = as.integer(),
    # Number of rounds (equivalent to days) 
    round = as.integer(),
    # Number of infectious people (diagnosed)
    new_inf_nominal = as.integer(), 
    # Cumulative infection cases (symptomatic + diagnosed)
    cic_nominal = as.numeric(),
    # Cumulative infection cases (not s)
    cic_total = as.numeric(),
    # Coef for the cox model
    coef = as.numeric(),
    # Se for the cox model
    se = as.numeric(),
    # Number of infected (not s) people in the placebo group
    placebo_inf = as.integer(),
    # Number of non-infected (s) people in the placebo group
    placebo_non_inf = as.integer(),
    # Number of infected (not s) people in the vaccine group
    vaccine_inf = as.integer(),
    # Number of non-infected (s) people in the vaccine group
    vaccine_non_inf = as.integer(),
    # 97.5% CI upper
    ci_upper = as.numeric(),
    # 2.5% CI lower
    ci_lower = as.numeric())
  
  # Read the ndata (Population data)
  setwd("../ndata")
  ndata <- open_dataset(
    sources = paste0("h", h, "_ndata.parquet"),
    format = "parquet") |> 
    collect() 
  
  # Make the new columns
  ndata <- ndata |>
    mutate(inf_exposure = 0,
           net_inf = 0)
  
  # Read the matrix 
  setwd("../mtx")
  xdata0 <- readMM(paste0("h", h, ".mtx"))
  
  # Update the session number
  ndata <- ndata |> mutate(session = h)
  
  # Id for initial E stage people
  initial_e_id <- sample(people_n, initial_e, replace = FALSE)
  
  # Update the state and state_end for initial E stage people
  ndata <- ndata |>
    mutate(state = if_else(id %in% initial_e_id, "e", state),
           state_end = if_else(id %in% initial_e_id, e_period, state_end))
  
  # Id for initial R stage people
  initial_r_id <- sample(
    setdiff(1:people_n, initial_e_id), initial_r, replace = FALSE)
  
  # Update the state and state_end for initial R stage people
  ndata <- ndata |>
    mutate(
      state = if_else(
        (id %in% initial_r_id) & (symptomatic == 1), "r_sym", state),
      state = if_else(
        (id %in% initial_r_id) & (symptomatic == 0), "r_asym", state))
  
  # Newly diagnosed cases (i_sym_post)
  new_inf_nominal = 0 
  # Cumulative infected cases (i_sym_post and r_sym)
  cic_nominal = 0
  # Cumulative infection cases (not s)
  cic_total = 0
  
  
  
}
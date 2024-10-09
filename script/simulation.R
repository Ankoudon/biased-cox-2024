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
per_trans <- 0.037
# E period
e_period <- 3
# R period 
r_period <- 0
# 0: Without permutation, 1: With permutation
tactic = 1
# Number of rounds (equivalent to days) 
round_num <- 210
# Number of participants 
subjects_n <- 300000
# Length of time for vaccine to be effective  
sample_period <- 42
# Vaccine efficacy 
# 0.677, 0.669, 0.740, 0.897, 0.916, 0.941, 0.950
v_benefit <- 0.941
# External infection per week
external_inf <- 700
# Initial E stage people 
initial_e <- 1100
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
  
  for (j in 1:round_num) {
    
    # Update round
    ndata <- ndata |>
      mutate(round = j)
    
    if (tactic == "1") {
      
      test_sample_indices <- ndata$test_sample == 1
      
      shuffled_data <- ndata |>
        filter(test_sample == 1) |>
        slice_sample(n = sum(test_sample_indices))
      
      ndata[test_sample_indices, ] <- shuffled_data
      
    }
    
    # Update the survival days
    ndata <- ndata |>
      mutate(survd = if_else(
        (test_sample == 1) & (state == "s"),
        survd + 1, survd))
    
    # Update the vaccine effect
    ndata <- ndata |>
      mutate(
        true_effect = if_else(test_sample == 1 & vaccine_status == 1,
                              v_benefit, true_effect))
    
    # Update "new_i"
    # {0: otherwise,
    #  1: stage change from E to I,
    #  2: stage change from I_Sym_Pre to I_Sym_Post}
    ndata <- ndata |> mutate(
      new_i = if_else((state == "e") & (round == state_end), 1,
                      new_i))
    # Update E to I_Asym 
    ndata <- ndata |> mutate(
      state = if_else((new_i == 1) & (symptomatic == 0),
                      "i_asym", state),
      state_end = if_else((state == "i_asym") & (new_i == 1),
                          (round + i_period), state_end))
    
    # Update E to I_Sym_Pre 
    ndata <- ndata |> mutate(
      state = if_else((new_i == 1) & (symptomatic == 1),
                      "i_sym_pre", state),
      state_end = if_else((state == "i_sym_pre") & (new_i == 1),
                          (round + pre_sym_period), state_end))
    
    # Update "new_i"
    ndata <- ndata |> mutate(
      new_i = if_else((state == "i_sym_pre") & (round == state_end), 2, new_i))
    
    # Update "I_Sym_Pre" to "I_Sym_Post"
    ndata <- ndata |> mutate(
      state = if_else(new_i == 2, "i_sym_post", state),
      state_end = if_else((state == "i_sym_post") & (new_i == 2),
                          (round + i_period - pre_sym_period), state_end))
    
    # Update "new_r"
    # {0: otherwise, 1: stage change from I to R}
    ndata <- ndata |> mutate(
      new_r = if_else((state %in% c("i_sym_post", "i_asym")) &
                        (round == state_end), 1, new_r))
    
    # Update I_Sym_Post to R_Sym
    ndata <- ndata |> mutate(
      state = if_else((new_r == 1) & (state == "i_sym_post"), "r_sym", state),
      state_end = if_else((state == "r_sym") & (new_r == 1),
                          r_period, state_end))
    
    # Update I_Asym to R_Asym
    ndata <- ndata |> mutate(
      state = if_else((new_r == 1) & (state == "i_asym"), "r_asym", state),
      state_end = if_else((state == "r_asym") & (new_r == 1),
                          r_period, state_end))
    
    # Extract "state"
    state_vec <- ndata |> select(state) |> pull()
    
    # Calculate the number of contacts (matrix)
    contacts_mtx <- xdata0 %*% if_else(
      state_vec %in% c("i_asym", "i_sym_pre", "i_sym_post"), 1, 0)
    
    # Number of contacts (vector)
    contacts_vec <- contacts_mtx[, 1]
    
    # Update the transmissibility
    ndata <- ndata |> mutate(
      contacts = contacts_vec,
      required_contacts = (
        rgeom(people_n, per_trans * (1 - true_effect)) + 1))
    
    # Update the "new_e"
    ndata <- ndata |> mutate(
      new_e = if_else(
        (state == "s") & (required_contacts <= contacts), 1, new_e))
    
    # Update the inf_exposure
    ndata <- ndata |>
      mutate(
        inf_exposure = if_else(
          (test_sample == 1) & (state == "s") & (new_e == 0),
          inf_exposure + contacts_vec,
          if_else( (test_sample == 1) & (state == "s") & (new_e == 1),
                   inf_exposure + required_contacts, inf_exposure)))
    
    # Update S to E
    ndata <- ndata |> mutate(
      state = if_else(new_e == 1, "e", state),
      state_end = if_else((state == "e") & (new_e == 1),
                          (round + e_period), state_end))
    
    # Save the exposure data from 143rd round
    if (j == (91 + sample_period)) {
      
      main_exposure_data <- ndata |> filter(test_sample == 1) |>
        select(id, round, vaccine_status, inf_exposure, degree, state)
      
    }
    
    if (j > (91 + sample_period)) {
      
      part_exposure_data <- ndata |> filter(test_sample == 1) |>
        select(id, round, vaccine_status, inf_exposure, degree, state)
      
      main_exposure_data <- add_row(main_exposure_data, part_exposure_data)
      
    }
    
    # External infection per week
    if (j %% 7 == 0) {
      
      # Extract candidates for external infection
      external_inf_id <- ndata |>
        filter(state == "s") |>
        sample_n(external_inf) |>
        pull(id)
      
      # Determine the external infections
      ndata <- ndata |>
        mutate(new_e = if_else((id %in% external_inf_id) &
                                 (true_effect < runif(people_n)), 1, new_e))
      
      # Update S to E
      ndata <- ndata |> 
        mutate( 
          state = if_else(new_e == 1, "e", state), 
          state_end = if_else((state == "e") & (new_e == 1),
                              (round + e_period), state_end))
      
    }
    
    # Update new_inf_nominal
    new_inf_nominal <- ndata |> select(new_i) |> filter(new_i == 2) |> nrow() 
    
    # Update cic_nominal
    cic_nominal <- ndata |> select(state) |>
      filter(state %in% c("i_sym_post", "r_sym")) |> nrow()
    
    # Update cic_total
    cic_total <- ndata |> select(state) |> filter(state != "s") |> nrow()
    
    # Initialize "new_e", "new_i", and "new_r"
    ndata <- ndata |> mutate(new_e = 0, new_i = 0, new_r = 0)
    
    # Conduct the second dose
    ndata <- ndata |> mutate(test_dose2 = if_else(
      (round == dose2_date) & (state %in% c("s")), 1, test_dose2))
    
    
    # Get the test sample
    ndata <- ndata |> mutate(
      test_sample = if_else(
        (round == sample_date) & (state == "s") & (test_dose2 == 1),
        1, test_sample))
    
    # Collect the participants for the clinical trial
    subjects = NULL
    
    if (j == 91) {
      
      # Choose the vaccine subjects
      subjects <- ndata |>
        filter(((state == "s")) & (vaccine_status == -1)) |>
        sample_n(subjects_n) |>
        pull(id)
      
      
    } 
    
    if (length(subjects) == subjects_n) {
      
      # Choose the placebo subjects
      placebo_med = sample(subjects, subjects_n/2, replace = FALSE)
      # Choose the vaccine subjects
      vaccine_med = setdiff(subjects, placebo_med)
      
      # Update for the subject group
      ndata <- ndata |>
        mutate(
          vaccine_status = if_else(
            (id %in% vaccine_med), 1,
            if_else((id %in% placebo_med), 0, vaccine_status)),
          vaccine_date = if_else(
            (id %in% subjects), round, vaccine_date),
          sample_date = if_else(
            (id %in% subjects), (round + sample_period), sample_date),
          dose2_date = if_else((id %in% subjects), (round + 28), dose2_date))
    }
    
    # Values for the result tibble
    cox_coef <- NA
    cox_se <- NA
    pla_inf <- NA
    pla_non <- NA
    vac_inf <- NA
    vac_non <- NA
    ci_upper <- NA
    ci_lower <- NA
    
    # Check if there are infected people in the sample
    sample_inf <- ndata |> filter((state != "s") & (test_sample == 1)) |> nrow()
    
    # Cox hazard model
    if (sample_inf > 0) {
      
      # Collect the test sample
      cox_data <- ndata |>
        filter(test_sample == 1) |>
        select(id, round, survd, vaccine_status, state)
      
      # Assign the censor flag
      cox_data <- cox_data |>
        mutate(censor = if_else(state == "s", 0, 1))

      cox_test <- coxphf(formula = Surv(
        time = survd, event = censor
      ) ~ vaccine_status, data = cox_data)
      
      # Coef for the cox
      cox_coef <- cox_test$coefficients[[1]]
      
      # SE for the coef
      cox_se <- sqrt(diag(cox_test$var))
      
      # Confidence interval
      ci_upper <- cox_test$ci.upper[[1]]
      ci_lower <- cox_test$ci.lower[[1]]
      
      # Number of infected people in the placebo group
      pla_inf <- cox_data |>
        filter((vaccine_status == 0) & (censor == 1)) |>
        nrow()
      
      # Number of non-infected people in the placebo group
      pla_non <- cox_data |>
        filter((vaccine_status == 0) & (censor == 0)) |>
        nrow()
      
      # Number of infected people in the vaccine group
      vac_inf <- cox_data |>
        filter((vaccine_status == 1) & (censor == 1)) |>
        nrow()
      
      # Number of non-infected people in the vaccine group
      vac_non <- cox_data |>
        filter((vaccine_status == 1) & (censor == 0)) |>
        nrow()
      
    }
    
    # Update "result"
    new_result <- tibble(
      session = h,
      round = j,
      new_inf_nominal = new_inf_nominal,
      cic_nominal = cic_nominal,
      cic_total = cic_total,
      coef = cox_coef,
      se = cox_se,
      placebo_inf = pla_inf,
      placebo_non_inf = pla_non,
      vaccine_inf = vac_inf,
      vaccine_non_inf = vac_non,
      ci_upper = ci_upper,
      ci_lower = ci_lower)
    
    result <- add_row(result, new_result)
    
    # Save the data
    if (j == round_num) {
      
      # Save the general result data
      setwd("../result")
      write_parquet(
        as_arrow_table(result),
        paste0("h", h, "_", v_benefit,"_result", "_t", tactic, ".parquet"))
      
      # Save the population data
      setwd("../end_ndata")
      write_parquet(
        as_arrow_table(ndata),
        paste0("h", h, "_", v_benefit, "_ndata", "_t", tactic, ".parquet"))
      
      # Save the exposure data
      setwd("../exposure")
      write_parquet(
        as_arrow_table(main_exposure_data),
        paste0("h", h, "_", v_benefit, "_exposure", "_t", tactic, ".parquet"))}}
  
  NULL
  
  }
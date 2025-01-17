################################################################################

# Thesis: Impact of Correlation of Unmeasured Exposure Status across Time on the
# Estimation of Vaccine Efficacy
# Periods On Estimating Vaccine Efficacy", Hiroyasu Ando, A. James O’Malley, Akihiro Nishi
# Purpose: Analyze the data

########## Table 1, VE_SAR, and VE_IR ##########

ve_sar_vec <- c()
ve_ir_vec <- c()

vac_non_inf <- 0
vac_inf <- 0
pla_non_inf <- 0
pla_inf <- 0

path_2 <- paste0("~/Desktop/biased-cox-2024/end_ndata/h")

# 0: Without permutation, 1: all-but-participants permutation, 2: all-but-bubbles permutation
tactic <- 0
# Vaccine efficacy 
v_benefit <- 0.941

for (h in 1:100) {
  
  print(h)

  data_end <- open_dataset(
    sources = paste0(
      path_2, h, "_", v_benefit, "_ndata", "_t", tactic, ".parquet"),
    format = "parquet") |> collect() 
  
  # Vaccine group
  vaccine <- data_end |>
    filter(vaccine_status == 1 & test_sample == 1)
  
  # Infected cases in the vaccine group
  vaccine_inf <- vaccine |>
    filter(state != "s")
  
  # Exposure in the vaccine group
  vaccine_expo <- vaccine |>
    pull(inf_exposure) 
  
  # Placebo group
  placebo <- data_end |>
    filter(vaccine_status == 0 & test_sample == 1)
  
  # Infected cases in the placebo group
  placebo_inf <- placebo |>
    filter(state != "s")
  
  # Exposure in the placebo group
  placebo_expo <- placebo |>
    pull(inf_exposure)
  
  # Count the number of non-infected cases in the vaccine group
  vac_non_inf <- vac_non_inf + (vaccine |> filter(state == "s") |> nrow())
  # Count the number of infected cases in the vaccine group
  vac_inf <- vac_inf + (vaccine_inf |> nrow())
  # Count the number of non-infected cases in the placebo group
  pla_non_inf <- pla_non_inf + (placebo |> filter(state == "s") |> nrow())
  # Count the number of infected cases in the placebo group
  pla_inf <- pla_inf + (placebo_inf |> nrow())
  
  # Calculate the vaccine efficacy (sar)
  ve_sar <- 1 - (
    (vaccine_inf |> nrow())/sum(vaccine_expo)) / 
    ((placebo_inf |> nrow())/sum(placebo_expo))
  # Calculate the vaccine efficacy (ir)
  ve_ir <- 1 - (
    (vaccine_inf |> nrow())/sum(vaccine$survd)) / 
    ((placebo_inf |> nrow())/sum(placebo$survd))
  
  # Append the vaccine efficacy to the vector
  ve_sar_vec <- c(ve_sar_vec, ve_sar)
  ve_ir_vec <- c(ve_ir_vec, ve_ir)
  
}

# Calculate the quantiles
quantile(ve_sar_vec, c(0.25, 0.5, 0.75))
quantile(ve_ir_vec, c(0.25, 0.5, 0.75))

# Create the table 1
matrix(c(vac_inf, pla_inf, vac_non_inf, pla_non_inf), nrow = 2)

########## Table 2 and Table 3 ##########

non_expo_non_inf <- 0
non_expo_inf <- 0
expo_non_inf <- 0
expo_inf <- 0

non_expo1_non_expo2 <- 0
non_expo1_expo2 <- 0
expo1_non_expo2 <- 0
expo1_expo2 <- 0

path_2 <- paste0("~/biased-cox-2024/exposure/h")

# 0: Without permutation, 1: all-but-participants permutation, 2: all-but-bubbles permutation
tactic <- 0
# Vaccine efficacy 
v_benefit <- 0.941

for (h in 1:100) {
  
  print(h)
  
  data <- open_dataset(
    sources = paste0(
      path_2, h, "_", v_benefit, "_exposure", "_t", tactic, ".parquet"),
    format = "parquet") |> 
    collect() 
  
  data <- data |> arrange(id, round)
  
  # Calculate the daily exposure
  data <- data |> 
    mutate(daily_exposure = inf_exposure - lag(inf_exposure)) |> 
    filter(((round >= 134) & (state == "s")) | 
             ((round >= 134) & (state == "e") & (daily_exposure > 0)))
  
  # Revise the daily exposure for the first day (134)
  data <- data |> 
    mutate(daily_exposure = ifelse(round == 134,
                                   inf_exposure,
                                   daily_exposure))
  
  # Count the number of non-infected cases with no exposure
  new_non_expo_non_inf <- data |> 
    filter(daily_exposure == 0 & state == "s") |> nrow()
  
  # Count the number of infected cases with no exposure
  new_non_expo_inf <- data |> 
    filter(state == "e" & daily_exposure == 0) |> nrow()
  
  # Count the number of non-infected cases with exposure
  new_expo_non_inf <- data |> 
    filter(daily_exposure > 0 & state == "s") |> nrow()
  
  # Count the number of infected cases with exposure
  new_expo_inf <- data |> 
    filter(state == "e" & daily_exposure > 0) |> nrow()
  
  # Update the counts
  non_expo_non_inf <- new_non_expo_non_inf + non_expo_non_inf
  non_expo_inf <- new_non_expo_inf + non_expo_inf
  expo_non_inf <- new_expo_non_inf + expo_non_inf
  expo_inf <- new_expo_inf + expo_inf
  
  # Calculate the daily exposure for the previous period
  data <- data |>
    mutate(daily_exposure_2 = lag(daily_exposure)) |> filter(round >= 135)
  
  # Count the number of cases which have no exposure in both periods
  new_non_expo1_non_expo2 <- data |> 
    filter(daily_exposure == 0 & daily_exposure_2 == 0) |> nrow()
  
  # Count the number of cases which have no exposure in this period and exposure in the previous period
  new_non_expo1_expo2 <- data |> 
    filter(daily_exposure == 0 & daily_exposure_2 > 0) |> nrow()
  
  # Count the number of cases which have exposure in this period and no exposure in the previous period
  new_expo1_non_expo2 <- data |>
    filter(daily_exposure > 0 & daily_exposure_2 == 0) |> nrow()
  
  # Count the number of cases which have exposure in both periods 
  new_expo1_expo2 <- data |>
    filter(daily_exposure > 0 & daily_exposure_2 > 0) |> nrow()
  
  # Update the counts
  non_expo1_non_expo2 <- new_non_expo1_non_expo2 + non_expo1_non_expo2
  non_expo1_expo2 <- new_non_expo1_expo2 + non_expo1_expo2
  expo1_non_expo2 <- new_expo1_non_expo2 + expo1_non_expo2
  expo1_expo2 <- new_expo1_expo2 + expo1_expo2
  
}

# Create the table 2
matrix(c(expo_inf, non_expo_inf, expo_non_inf, non_expo_non_inf), nrow = 2)

# Create the table 3
matrix(c(expo1_expo2, expo1_non_expo2, non_expo1_expo2, non_expo1_non_expo2),
       nrow = 2)

########## Fisher's combined probability test ##########

library(metap)

# p-values from 2x2 tables
p_vec_fisher <- c()

path_2 <- paste0("~/biased-cox-2024/exposure/h")

# 0: Without permutation, 1: all-but-participants permutation, 2: all-but-bubbles permutation
tactic <- 1
# Vaccine efficacy 
v_benefit <- 0.941

for (h in 1:100) {
  
  print(h)
  
  data <- open_dataset(
    sources = paste0(
      path_2, h, "_", v_benefit, "_exposure", "_t", tactic, ".parquet"),
    format = "parquet") |> 
    collect() 
  
  data <- data |> arrange(id, round)
  
  # Calculate the daily exposure
  data <- data |> 
    mutate(daily_exposure = inf_exposure - lag(inf_exposure)) |> 
    filter(((round >= 134) & (state == "s")) | 
             ((round >= 134) & (state == "e") & (daily_exposure > 0)))
  
  # Revise the daily exposure for the first day (134)
  data <- data |> 
    mutate(daily_exposure = ifelse(round == 134,
                                   inf_exposure,
                                   daily_exposure))
  
  # Calculate the daily exposure for the previous period
  data <- data |>
    mutate(daily_exposure_2 = lag(daily_exposure)) |>
    filter(round >= 135)
  
  # Create a 2x2 table every day
  for (i in 135:210) {
    
    print(i)
    
    data_day <- data |> filter(round == i)
    
    # Count the number of cases which have no exposure in both periods
    non_expo1_non_expo2 <- data_day |> 
      filter(daily_exposure == 0 & daily_exposure_2 == 0) |> nrow()
    
    # Count the number of cases which have no exposure in this period and exposure in the previous period
    non_expo1_expo2 <- data_day |> 
      filter(daily_exposure == 0 & daily_exposure_2 > 0) |> nrow()
    
    # Count the number of cases which have exposure in this period and no exposure in the previous period
    expo1_non_expo2 <- data_day |>
      filter(daily_exposure > 0 & daily_exposure_2 == 0) |> nrow()
    
    # Count the number of cases which have exposure in both periods
    expo1_expo2 <- data_day |>
      filter(daily_exposure > 0 & daily_exposure_2 > 0) |> nrow()
    
    # Make 2x2 table
    mtx <- matrix(c(expo1_expo2,
                    expo1_non_expo2,
                    non_expo1_expo2,
                    non_expo1_non_expo2), nrow = 2)
    # Fisher's exact test
    fisher_test <- fisher.test(mtx)
    # Append the p-value to the vector
    p_vec_fisher <- c(fisher_test$p.value, p_vec_fisher)
    
  }
  
}

# Fisher's combined probability test
sumlog(p_vec_fisher)



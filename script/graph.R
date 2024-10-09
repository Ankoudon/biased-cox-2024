################################################################################

# Thesis: "Impact of Correlation of Exposure between Two Consecutive
# Periods On Estimating Vaccine Efficacy", Hiroyasu Ando, Akihiro Nishi
# Purpose: Plot the estimated vaccine efficacy 

library(tidyverse)
library(dplyr)
library(arrow)
library(ggplot2)
library(gridExtra)

########## Graph A ##########

path_2 <- paste0("~/Desktop/biased-cox-2024/result/h")

# Vaccine efficacy 
v_benefit <- 0.941

# Read the first result data
# 0: Without permutation, 1: With permutation

data_0 <- open_dataset(
  sources = paste0(
    path_2, 1, "_", v_benefit, "_result", "_t", 0, ".parquet"),
  format = "parquet") |> collect() 

data_1 <- open_dataset(
  sources = paste0(
    path_2, 1, "_", v_benefit, "_result", "_t", 1, ".parquet"),
  format = "parquet") |> collect() 

# Put all the result data together
for (h in 2:100) {
  
  print(h)
  
  new_data_0 <- open_dataset(
    sources = paste0(
      path_2, h, "_", v_benefit, "_result", "_t", 0, ".parquet"),
    format = "parquet") |> collect() 
  
  data_0 <- add_row(data_0, new_data_0)
  
  new_data_1 <- open_dataset(
    sources = paste0(
      path_2, 1, "_", v_benefit, "_result", "_t", 1, ".parquet"),
    format = "parquet") |> collect() 
  
  data_1 <- add_row(data_1, new_data_1)
  
}

# Tables for the VE estimates each day
coef_table_0 <- tibble(
  day = as.integer(),
  coef = as.numeric())

coef_table_1 <- tibble(
  day = as.integer(),
  coef = as.numeric())


for (i in c(150, 165, 180, 195, 210)) {
  
  # Extract the VE estimates for each day
  new_coef_table_0 <- data_0 |> 
    filter(round == i) |> 
    mutate(coef = 1 - exp(coef),
           day = round) |>
    select(coef, day) 
  
  new_coef_table_1 <- data_1 |> 
    filter(round == i) |> 
    mutate(coef = 1 - exp(coef),
           day = round) |>
    select(coef, day)
  
  # Put the VE estimate tables together
  coef_table_0 <- add_row(coef_table_0, new_coef_table_0)
  coef_table_1 <- add_row(coef_table_1, new_coef_table_1)
  
}

coef_table_0 <- coef_table_0 |>
  mutate(type = "Without Node Permutation") 

coef_table_1 <- coef_table_1 |>
  mutate(type = "With Node Permutation") 


# Combine the tables
coef_table <- add_row(coef_table_0, coef_table_1) |>
  mutate(type = factor(
    type, levels = c("Without Node Permutation", "With Node Permutation")))

# Depict a graph
coef_graph <- coef_table |>
  ggplot(aes(x = factor(day), y = coef, fill = type)) + 
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.5) +  
  scale_fill_manual(
    values = c("With Node Permutation" = "#0072B2",
               "Without Node Permutation" = "#E69F00"), name = "") +
  geom_hline(yintercept = 0.941, color = "red", size = 1) +
  labs(x = "Day from the Beginning of Network-based Simulations",
       y = "Vaccine Efficacy (VE) Estimate") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.position = c(1.05, 0.01),  
        legend.justification = c(1, 0))

# Quantile of VE estimates (day 210)
quantile(coef_table |>
           filter(day == 210 & type == "Without Permutation") |>
           pull(coef), c(0.25, 0.5, 0.75))
quantile(coef_table |>
           filter(day == 210 & type == "With Permutation") |>
           pull(coef), c(0.25, 0.5, 0.75))

########## Graph B (Robustness) ##########

# Table for the VE estimates each true VE
ve_tibble <- tibble(VE = as.numeric(),
                    coef = as.numeric())

path_2 <- paste0("~/Desktop/biased-cox-2024/result/h")

for (i in c("0.677", "0.669", "0.740", "0.897", "0.916", "0.941", "0.950")) {
  
  # Read the first result data
  data <- open_dataset(
    sources = paste0(path_2, 1, "_", i, "_result", "_t", t, ".parquet"),
    format = "parquet") |> 
    collect() 
  
  # Put all the result data together
  for (h in 2:100) {
    
    print(h)
    
    new_data <- open_dataset(
      sources = paste0(path_2, h, "_", i, "_result", "_t", t, ".parquet"),
      format = "parquet") |> 
      collect() 
    
    data <- add_row(data, new_data)
    
  }
  
  data <- data |> filter(round == 210) |>
    mutate(VE = as.numeric(i),
           coef = ((1 - exp(coef)) - as.numeric(i))) |>
    select(VE, coef) 
  
  ve_tibble <- add_row(ve_tibble, data)
  
}


robust_graph <- ve_tibble |>
  ggplot(aes(x = factor(VE), y = coef, fill = factor(VE))) + 
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.5) +
  scale_fill_manual(values = c("0.669" = "#e41a1c",
                               "0.677" = "#377eb8",
                               "0.74" = "#4daf4a",
                               "0.897" = "#984ea3",
                               "0.916" = "#ff7f00",
                               "0.941" = "#ffff33",
                               "0.95" = "#a65628"), name = "",
                    labels = c("0.669" = "Ad26.COV2.S (0.669)",
                               "0.677" = "Coronavac (0.677)",
                               "0.74" = "ChAdOx1 nCoV-2019 (0.740)",
                               "0.897" = "NVX-CoV2373 (0.897)",
                               "0.916" = "Gam-COVID-Vac Sputnik V (0.916)",
                               "0.941" = "mRNA-1273 (0.941)",
                               "0.95" = "BNT162b2 (0.950)")) + 
  labs(x = "True VE", y = "VE Estimate - True VE") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.position = c(1.02, 0.01), 
        legend.justification = c(1, 0)) +
  scale_x_discrete(labels = c("0.669" = "0.669",
                              "0.677" = "0.677",
                              "0.74" = "0.740",
                              "0.897" = "0.897",
                              "0.916" = "0.916",
                              "0.941" = "0.941",
                              "0.95" = "0.950"))

########## Combine graph A and B ##########

ggarrange(coef_graph,
          robust_graph,
          labels = c("A", "B"),
          ncol = 2,
          nrow = 1,
          font.label = list(size = 16),
          label.x = c(0, 0))

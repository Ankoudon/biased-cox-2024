################################################################################

# Thesis: Impact of Correlation of Unmeasured Exposure Status across Time on the
# Estimation of Vaccine Efficacy
# Periods On Estimating Vaccine Efficacy", Hiroyasu Ando, A. James O’Malley, Akihiro Nishi
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
# 0: Without permutation, 1: all-but-participants permutation, 2: all-but-bubbles permutation

data_0 <- open_dataset(
  sources = paste0(
    path_2, 1, "_", v_benefit, "_result", "_t", 0, ".parquet"),
  format = "parquet") |> collect() 

data_1 <- open_dataset(
  sources = paste0(
    path_2, 1, "_", v_benefit, "_result", "_t", 1, ".parquet"),
  format = "parquet") |> collect() 

data_2 <- open_dataset(
  sources = paste0(
    path_2, 1, "_", v_benefit, "_result", "_t", 2, ".parquet"),
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
      path_2, h, "_", v_benefit, "_result", "_t", 1, ".parquet"),
    format = "parquet") |> collect() 
  
  data_1 <- add_row(data_1, new_data_1)
  
  new_data_2 <- open_dataset(
    sources = paste0(
      path_2, h, "_", v_benefit, "_result", "_t", 2, ".parquet"),
    format = "parquet") |> collect() 
  
  data_2 <- add_row(data_2, new_data_2)
  
}

# Tables for the VE estimates each day
coef_table_0 <- tibble(
  day = as.integer(),
  coef = as.numeric())

coef_table_1 <- tibble(
  day = as.integer(),
  coef = as.numeric())

coef_table_2 <- tibble(
  day = as.integer(),
  coef = as.numeric())


for (i in c(170, 190, 210)) {
  
  # Extract the VE estimates for each day
  new_coef_table_0 <- data_0 |> 
    filter(round == i) |> 
    mutate(coef = 1 - exp(coef),
           day = round) |>
    dplyr::select(coef, day) 
  
  new_coef_table_1 <- data_1 |> 
    filter(round == i) |> 
    mutate(coef = 1 - exp(coef),
           day = round) |>
    dplyr::select(coef, day)
  
  new_coef_table_2 <- data_2 |> 
    filter(round == i) |> 
    mutate(coef = 1 - exp(coef),
           day = round) |>
    dplyr::select(coef, day)
  
  # Put the VE estimate tables together
  coef_table_0 <- add_row(coef_table_0, new_coef_table_0)
  coef_table_1 <- add_row(coef_table_1, new_coef_table_1)
  coef_table_2 <- add_row(coef_table_2, new_coef_table_2)
  
}

coef_table_0 <- coef_table_0 |>
  mutate(type = "without permutation") 

coef_table_1 <- coef_table_1 |>
  mutate(type = "all-but-participants permutation") 

coef_table_2 <- coef_table_2 |>
  mutate(type = "all-but-bubbles permutation") 


# Combine the tables
coef_table <- bind_rows(coef_table_0, coef_table_1, coef_table_2) |>
  mutate(
    type = factor(
    type,
    levels = c("without permutation",
               "all-but-bubbles permutation",
               "all-but-participants permutation")))

# Depict a graph
coef_graph <- coef_table |>
  ggplot(aes(x = factor(day), y = coef, fill = type)) + 
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.5) +  
  scale_fill_manual(
    values = c("without permutation" = "#66B2FF",
               "all-but-bubbles permutation" = "#E69F00",
               "all-but-participants permutation" = "#FF69B4"), name = "") +
  #geom_hline(yintercept = 0, color = "red", size = 0.5) +
  geom_hline(yintercept = 0.941, color = "red", size = 0.5) +
  scale_y_continuous(limits = c(NA, 0.98),
                     breaks = seq(0.92, 0.96, by = 0.01)) +
  labs(x = "Days from the beginning of network-based simulations",
       y = "Vaccine efficacy (VE) estimate") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.position = c(1, 0.001),  
        legend.justification = c(1, 0))

# Quantile of VE estimates (day 210)
quantile(coef_table |>
           filter(day == 210 & type == "without permutation") |>
           pull(coef), c(0.25, 0.5, 0.75))
quantile(coef_table |>
           filter(day == 210 & type == "all-but-bubbles permutation") |>
           pull(coef), c(0.25, 0.5, 0.75))
quantile(coef_table |>
           filter(day == 210 & type == "all-but-participants permutation") |>
           pull(coef), c(0.25, 0.5, 0.75))


########## Graph B (Robustness) ##########

# Table for the VE estimates each true VE
ve_tibble <- tibble(VE = as.numeric(),
                    coef = as.numeric())

path_2 <- paste0("~/Desktop/biased-cox-2024/result/h")

for (i in c("0.677", "0.669", "0.74", "0.897", "0.916", "0.941", "0.95")) {
  
  # Read the first result data
  data <- open_dataset(
    sources = paste0(path_2, 1, "_", i, "_result", "_t", 0, ".parquet"),
    format = "parquet") |> 
    collect() 
  
  # Put all the result data together
  for (h in 2:100) {
    
    print(h)
    
    new_data <- open_dataset(
      sources = paste0(path_2, h, "_", i, "_result", "_t", 0, ".parquet"),
      format = "parquet") |> 
      collect() 
    
    data <- add_row(data, new_data)
    
  }
  
  data <- data |> filter(round == 210) |>
    mutate(VE = as.numeric(i),
           coef = ((1 - exp(coef)) - as.numeric(i))) |>
    dplyr::select(VE, coef) 
  
  ve_tibble <- add_row(ve_tibble, data)
  
}


robust_graph <- ve_tibble |>
  ggplot(aes(x = factor(VE), y = coef, fill = factor(VE))) + 
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.5) +
  scale_fill_manual(values = c("0.669" = "#001F3F",
                               "0.677" = "#003F7F",
                               "0.74" = "#005FBF",
                               "0.897" = "#007FFF",
                               "0.916" = "#3399FF",
                               "0.941" = "#66B2FF",
                               "0.95" = "#99CCFF"), name = "",
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
                              "0.95" = "0.950")) +
  geom_hline(yintercept = 0, color = "red", size = 0.5) 

########## Combine graph A and B ##########

ggarrange(coef_graph,
          robust_graph,
          labels = c("A", "B"),
          ncol = 2,
          nrow = 1,
          font.label = list(size = 16),
          label.x = c(0, 0))

########## Appendix ##########

path_2 <- paste0("~/Desktop/biased-cox-2024/result/h")

# Data from The New York Times, based on reports from state and local health agencies.
us_data <- as_tibble(read_csv("~/Desktop/biased-cox-2024/real_data/Country.csv"))

us_data <- us_data |>
  filter(as.Date('2020-05-01') <= date & date <= as.Date('2020-12-1'))
# 331,449,281 (Population, Census, April 1, 2020)
us_data <- us_data |> mutate(cases = cases / 331449281)
start_date <- as.Date('2020-04-30')

us_data <- us_data |>
  mutate(
    lower_bound = cases * 0.7,
    upper_bound = cases * 1.3
  )

gg <- ggplot() +
  labs(
    x = "Date",
    y = "Cumulative Incidence",
    color = "Legend") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

gg <- gg +
  geom_ribbon(data = us_data,
              aes(x = date, ymin = lower_bound, ymax = upper_bound), 
              fill = "#FFB3B3", alpha = 0.3) +
  geom_line(data = us_data, aes(x = date, y = cases, color = "Real Data ± 30%"), 
            linewidth = 1.5)

for (h in 1:100) {
  
  rdata <- open_dataset(
    sources = paste0(path_2, h, "_", v_benefit, "_result", "_t", 0, ".parquet"),
    format = "parquet") |> 
    collect() |>
    mutate(cases = cic_nominal / 10000000) |>
    filter(round <= 210)
  
  rdata <- rdata |>
    mutate(date = start_date + round)
  
  gg <- gg + geom_line(
    data = rdata, aes(x = date, y = cases, color = "Network-based Simulations"),
    linewidth = 0.1, alpha = 0.3)

}


gg <- gg + 
  scale_color_manual(values = c("Real Data ± 30%" = "#E74C3C",
                                "Network-based Simulations" = "#3498DB")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%Y-%m") +
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.position = c(0.2, 0.8),  
        legend.justification = c(0.5, 0.5),
        axis.title.x = element_text(size = 18),  
        axis.title.y = element_text(size = 18),  
        axis.text = element_text(size = 15),     
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.text = element_text(size = 15),  
        legend.title = element_blank())





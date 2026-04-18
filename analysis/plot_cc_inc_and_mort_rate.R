#* Title: plot_cc_inc_&_mort_rate
#* 
#* Code function: 
#*    Plot the incidence and mortality rate in Mexico from 2000 to 2010.
#*    
#*    The data was digitized from Figure 2 of the following source:
#*    https://pdfs.semanticscholar.org/785b/2991f2205aa84b3a0e580fa7ae7d4a64323d.pdf
#* 
#* Creation date: February 23 2026
#* Author: David U. Garibay-Treviño

# 01 Initial Setup --------------------------------------------------------

## 01.01 Clean environment ------------------------------------------------
remove(list = ls())

#* Refresh environment memory
gc()

## 01.02 Load libraries ----------------------------------------------------
library(ggplot2)
library(dplyr)



# 02 Load data ------------------------------------------------------------

df_inc_rate  <- read.csv(file = "data-raw/age-adjusted-inc-rate-cc.csv")
df_mort_rate <- read.csv(file = "data-raw/age-adjusted-mort-rate-cc.csv")



# 03 Merge data -----------------------------------------------------------

df_cc_rates <- bind_rows(df_inc_rate,
                         df_mort_rate,
                         .id = "Type")


# 03 Plotting -------------------------------------------------------------

n_base_size <- 18

plt_inc_rate <- ggplot(data = df_inc_rate,
       mapping = aes(x = year,
                     y = rate,
                     color = age_group,
                     group = age_group)) + 
  theme_bw(base_size = n_base_size) + 
  geom_line() + 
  geom_point() +
  scale_x_continuous(breaks = seq(1990, 2050, 2)) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  labs(x = "Year",
       y = "Mortality rate (per 100,000)",
       color = "Age group") +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 70))


plt_mort_rate <- ggplot(data = df_mort_rate,
       mapping = aes(x = year,
                     y = rate,
                     color = age_group,
                     group = age_group)) + 
  theme_bw(base_size = n_base_size) + 
  geom_line() + 
  geom_point() +
  scale_x_continuous(breaks = seq(1990, 2050, 2)) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  labs(x = "Year",
       y = "Mortality rate (per 100,000)",
       color = "Age group") +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 60))

plt_cc_rates <- ggplot(data = df_cc_rates,
       mapping = aes(x = year,
                     y = rate,
                     color = age_group,
                     group = interaction(age_group, type))) + 
  theme_bw(base_size = n_base_size) + 
  geom_line() + 
  geom_point() +
  scale_x_continuous(breaks = seq(1990, 2050, 2)) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  labs(x = "Year",
       y = "Rate per 100,000",
       color = "Age group") +
  coord_cartesian(ylim = c(0, 60)) +
  theme(legend.position = "bottom") +
  facet_wrap(~type)


# 04 Save plots -----------------------------------------------------------

n_width <- 10
n_height <- 6

ggsave(plot = plt_inc_rate,  filename = "figs/plt_inc_rate.png",  width = n_width, height = n_height)
ggsave(plot = plt_mort_rate, filename = "figs/plt_mort_rate.png", width = n_width, height = n_height)
ggsave(plot = plt_cc_rates,  filename = "figs/plt_cc_rates.png",  width = n_width*1.5, height = n_height)

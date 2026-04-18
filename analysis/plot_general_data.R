#* Title: plot_general_data
#* 
#* Code function:
#*    This code generates the plots for the following data
#*    - HPV prevalence in Mexico (https://www.sciencedirect.com/science/article/pii/S2667193X25001668?via%3Dihub)
#*    - 
#* 
#* 
#* Creation date: February 25 2026
#* Author: author_name

# 01 Initial Setup --------------------------------------------------------

## 01.01 Clean environment ------------------------------------------------
remove(list = ls())

#* Refresh environment memory
gc()

## 01.02 Load libraries ----------------------------------------------------
library(dplyr)
library(ggplot2)
library(readr)



# 02 Read data ------------------------------------------------------------

#* *HPV prevalence in Mexico*
#* source: Appendix of https://www.sciencedirect.com/science/article/pii/S2667193X25001668?via%3Dihub
df_HPV_prev_00 <- read_csv(file = "data-raw/HPV_overall_prev_mx_2023.csv")

#* *Age-standardized CC incidence rate*
#* source: digitized Fig 2 of https://pdfs.semanticscholar.org/785b/2991f2205aa84b3a0e580fa7ae7d4a64323d.pdf
df_cc_inc_rate <- read_csv(file = "data-raw/age-adjusted-inc-rate-cc.csv")

#* *Age-standardized CC mortality rate*
#* source: digitized Fig 2 of https://pdfs.semanticscholar.org/785b/2991f2205aa84b3a0e580fa7ae7d4a64323d.pdf
df_cc_mort_rate <- read_csv(file = "data-raw/age-adjusted-mort-rate-cc.csv")



# 03 Wrangle data ---------------------------------------------------------


## 03.01 HPV prevalence ----------------------------------------------------

# Add ordered levels to the age groups
df_HPV_prev_00 <- df_HPV_prev_00 %>% 
  mutate(`Age group` = factor(`Age group`,
                              levels = c("<25", "25-34", "35-44", "45-54",
                                         "55-64","65-74", "≥75"),
                              ordered = T))


# 03.02 CC incidence rate -------------------------------------------------

# Add ordered levels to age groups and years
df_cc_inc_rate <- df_cc_inc_rate %>% 
  mutate(age_group = factor(x = age_group,
                            levels = c("25-44", "45-49", "50-59", "60-64", "65+"),
                            ordered = T),
         year = as.character(year)) %>% 
  mutate(year = factor(year,
                       levels = seq(2000, 2010, 2),
                       ordered = T))

# 03.03 CC mortality rate -------------------------------------------------

# Add ordered levels to age groups and years
df_cc_mort_rate <- df_cc_mort_rate %>% 
  mutate(age_group = factor(x = age_group,
                            levels = c("25-44", "45-49", "50-59", "60-64", "65+"),
                            ordered = T),
         year = as.character(year)) %>% 
  mutate(year = factor(year,
                       levels = seq(2000, 2010, 2),
                       ordered = T))

# 04 Plotting -------------------------------------------------------------

n_base_size <- 18


plt_HPV_prev <- ggplot(data = df_HPV_prev_00,
                       mapping = aes(x = `Age group`,
                                     y = prevalence,
                                     ymin = `Lower CI`,
                                     ymax = `Upper CI`)) + 
  theme_bw(base_size = n_base_size) +
  geom_point() +
  geom_errorbar() + 
  labs(y = "Overall HPV prevalence")
  

plt_CC_inc <- ggplot(data = df_cc_inc_rate,
                       mapping = aes(x = year,
                                     y = inc_rate,
                                     group = age_group,
                                     color = age_group,
                                     shape = age_group)) + 
  theme_bw(base_size = n_base_size) +
  geom_point() + 
  geom_line() + 
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Year",
       y = "Age-standardized incidence rate\nper 100,000",
       color = "Age group",
       shape = "Age group") +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom")


plt_CC_mort <- ggplot(data = df_cc_mort_rate,
                       mapping = aes(x = year,
                                     y = rate,
                                     group = age_group,
                                     color = age_group,
                                     shape = age_group)) + 
  theme_bw(base_size = n_base_size) +
  geom_point() + 
  geom_line() + 
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Year",
       y = "Age-standardized mortality rate\nper 100,000",
       color = "Age group",
       shape = "Age group") +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom")


# 05 Save outputs ---------------------------------------------------------


n_width <- 10
n_height <- 6

ggsave(plot = plt_HPV_prev, filename = "figs/plt_HPV_prev.png", width = n_width, height = n_height)
ggsave(plot = plt_CC_inc,   filename = "figs/plt_CC_inc.png",   width = n_width, height = n_height)
ggsave(plot = plt_CC_mort,  filename = "figs/plt_CC_mort.png",  width = n_width, height = n_height)

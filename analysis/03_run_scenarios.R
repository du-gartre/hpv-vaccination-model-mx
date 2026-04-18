#* Title: 03_run_scenarios
#* 
#* Code function: This code uses the twig model and the calibrated parameters
#*    from `analysis/02_calibration.R` to run different vaccination 
#*    scenarios:
#*    
#*    1. No Treatment
#*    2. Vaccinate only women
#*       - 40% coverage
#*       - 60% coverage
#*       - 80% coverage 
#*    3. Gender neutral vaccination
#*       - 40% coverage
#*       - 60% coverage
#*       - 80% coverage
#* 
#* Creation date: April 14 2026
#* Author: David U. Garibay Treviño

# 01 Initial Setup --------------------------------------------------------

## 01.01 Clean environment ------------------------------------------------
remove(list = ls())

#* Refresh environment memory
gc()

## 01.02 Load libraries ----------------------------------------------------
library(twig)
library(readr)
library(dplyr)
library(readxl)
library(patchwork)

source(file = "R/twig_funs.R")

# 02 Load data ------------------------------------------------------------

## 02.01 Mortality rates ---------------------------------------------------

#* Load the mortality rates/probs for mexican females, during 2020
#* Mortality rates go from 0 to 109 years of age
df_life_tables <- read_csv(file = "data/df_life_table.csv") %>% 
  # Keep data for females older than 16 years of age
  filter(age >= 16)

# Extract vector with mortality rates by age
v_r_mort <- df_life_tables$mort_rate

names(v_r_mort) <- paste0("age_", 16:109)

n_cycles <- length(v_r_mort)




# 03 Model parameters -----------------------------------------------------



## 03.01 General parameters ------------------------------------------------

# Number of simulations
n_sim <- 1e3


# Define model parameters using calibrated values
df_params <- dplyr::tibble(
  
  # # Proportion of vaccinated people by age
  # prop_vac_women = 0.8, #* *Assumption*
  # prop_vac_men   = 0.4, #* *Assumption*
  
  # Transition parameters parameters
  r_H_II     = 0.016260416,  #* *calibrated value* check script `analysis/02.01_calibration.R`
  r_II_PI    = 0.005949991,  #* *calibrated value* check script `analysis/02.01_calibration.R`
  r_PI_PL    = 0.010190967,  #* *calibrated value* check script `analysis/02.01_calibration.R`
  r_PL_CL    = 0.006369672,  #* *calibrated value* check script `analysis/02.01_calibration.R`
  cc_c_mort_rate = 6.58/1e5, # CC crude mortality rate. Source: section 3.3.3 of Human Papillomavirus and Related Diseases Report MEXICO
  hr_PL      = 1.1,  #* *assumed values* HR applied to baseline mort rate when in PL
  hr_CL      = 1.3,  #* *assumed values* HR applied to baseline mort rate when in CL
  hr_hpv_vac = 0.05, #* HR representing protective effect of vaccination among females (assuming 95% of efficacy, Baussano, et al., 2017)
  hr_hpv_men = 0.2,  #* HR representing protective HPV herd effect of gender neutral vaccination
  
  # Recovery rates
  r_II_HV = 0.7, # source: https://pmc.ncbi.nlm.nih.gov/articles/PMC6223532/#:~:text=The%20lrHPVs%20are%20mainly%20associated,20%25%20of%20cervical%20cancers%20worldwide.&text=Approximately%2070%25%20of%20HPV%20infections,in%20the%20remainder%20of%20cases.&text=Although%20persistent%20infection%20with%20hrHPV,effective%20cell%2Dmediated%20immune%20response.&text=Therefore%2C%20HIV%2Dpositive%20individuals%20who,risk%20factors%20among%20Nigerian%20women.
  r_PI_HV = 0.5, #* *check* source: https://link.springer.com/article/10.1186/s12905-023-02764-8
  r_PL_HV = 0.4, #* *assumed values* (ordinal recovery rate)
  r_CL_HV = 0.2, #* *assumed values* (ordinal recovery rate)
  
  # Cost parameters (USD)
  #* https://www.saludpublica.mx/index.php/spm/article/view/7374
  #* https://pmc.ncbi.nlm.nih.gov/articles/PMC3025113/#Sec7
  #* https://www.scielosp.org/pdf/spm/2014.v56n5/502-510#:~:text=Objective.,of%20undetected%20cervical%20cancer%20cases.&text=cervical%20cancer;%20Mexico-,Granados%2DGarc%C3%ADa%20V%2C%20Flores%20YN%2C%20P%C3%A9rez%20R%2C%20Rudolph,;56:502%2D510.
  c_vacc = ((331*2)*1.13)/17.95, # Using cost of CENETEC 2023, considering inflation to Jan 2026 and MXN-USD exchange rate of Jan 2026
  # c_trt_FOV  = n_women*prop_vac_women*c_vacc,
  # c_trt_vacc = c_vacc*(n_women*prop_vac_women + n_men*prop_vac_men),
  
  c_HV    = c_vacc + 19.5,     #* Vaccination + annual pap smear (Beal et al, 2014). June 2013 price converted to Jan 2026 assuming general CPI inflation
  c_HU    = 19.5,              #* Assuming annual pap smear (Beal et al, 2014). June 2013 price converted to Jan 2026 assuming general CPI inflation
  c_II   = (2000*2.52)/17.95,  #* Cost of genital wart treatment (Insinga, et al, 2007) MXN price from Jan 2005 to Jan 2026, then applied official MXN-USD exchange rate of Jan 2026
  c_PI   = (2000*2.52)/17.95,  #* Cost of genital wart treatment (Insinga, et al, 2007) MXN price from Jan 2005 to Jan 2026, then applied official MXN-USD exchange rate of Jan 2026
  c_PL   = (15849*2.52)/17.95, #* Cost of CIN treatment (Insinga, et al, 2007) MXN price from Jan 2005 to Jan 2026, then applied official MXN-USD exchange rate of Jan 2026
  # c_CL   = 952.54*12,        #* Mean annual cost of CC treatment - all stages (Granados-García et al., 2019). June 2017 price converted to Jan 2026 assuming general CPI inflation
  c_CL   = 11000,              #* Mean annual cost of CC treatment - all stages (Granados-García et al., 2019). June 2017 price converted to Jan 2026 assuming general CPI inflation
  c_D    = 0,                  #* Assuming the public healthcare provider
  
  # # Utility parameters
  u_HV = 1,               # Time-dependent utility weight for healthy state 
  u_HU = 1,               # Time-dependent utility weight for healthy state
  u_II = 0.91,            # Utility weight for genital wart (Elbasha et al. 2007)
  u_PI = 0.91,            # Utility weight for CIN 1 (Elbasha et al. 2007)
  u_PL = 0.87,            # Utility weight for CIN 2/3 (Elbasha et al. 2007)
  u_CL = (0.76 + 0.67)/2  # Average Utility weight for localized and regional cervical cancer (Elbasha et al. 2007)
)


## 03.02 Scenario-specific parameters --------------------------------------


# assuming 40% coverage and 95% of efficacy
df_params_40 <- df_params
df_params_40$hr_hpv_men <- 1/(0.67/(0.4*0.95)) # Value calculated using values from Table 2 of (Baussano, et al., 2017)

# assuming 60% coverage and 95% of efficacy
df_params_60 <- df_params
df_params_60$hr_hpv_men <- 1/(0.90/(0.6*0.95)) # Value calculated using values from Table 2 of (Baussano, et al., 2017)

# assuming 80% coverage and 95% of efficacy
df_params_80 <- df_params
df_params_80$hr_hpv_men <- 1/(1.00/(0.8*0.95)) # Value calculated using values from Table 2 of (Baussano, et al., 2017)


# 04 Define twigs ---------------------------------------------------------

#* Define twigs for different starting scenarios:

#* 1. Twig for 40% vaccination coverage
hpv_twig_40 <- twig() + 
  decisions(names = c(no_trt, FOV, GNV)) + 
  states(names = c(HV, HU, HPV_II, HPV_PI, PL, CL, Death),
         init_probs = c(0.4, 0.6, 0, 0, 0 , 0, 0),
         max_cycles = c(1, 1, 1, 1, 1, 1, 1)) +
  event(name = death_event,
        options = c(yes, none),
        probs = c(p_die, leftover), 
        transitions = c(Death, second_event)) +
  event(name = second_event,
        options     = c(recover,   get_II,   get_PI,   get_PL,   get_CL,   none),
        probs       = c(p_recover, p_get_II, p_get_PI, p_get_PL, p_get_CL, leftover),
        transitions = c(HV,        HPV_II,   HPV_PI,   PL,       CL,       stay)) +
  payoffs(names = c(cost, utility, incidence),
          #* Discount rate assuming annual cycle length
          #* use 5% for LMIC from (Haacker, Hallett & Atun, 2020)
          discount_rates = c(0.05, 0.05, 0))

#* 2. Twig for 60% vaccination coverage
hpv_twig_60 <- twig() + 
  decisions(names = c(no_trt, FOV, GNV)) + 
  states(names = c(HV, HU, HPV_II, HPV_PI, PL, CL, Death),
         init_probs = c(0.6, 0.4, 0, 0, 0 , 0, 0),
         max_cycles = c(1, 1, 1, 1, 1, 1, 1)) + 
  event(name = death_event,
        options = c(yes, none),
        probs = c(p_die, leftover), 
        transitions = c(Death, second_event)) +
  event(name = second_event,
        options     = c(recover,   get_II,   get_PI,   get_PL,   get_CL,   none),
        probs       = c(p_recover, p_get_II, p_get_PI, p_get_PL, p_get_CL, leftover),
        transitions = c(HV,        HPV_II,   HPV_PI,   PL,       CL,       stay)) +
  payoffs(names = c(cost, utility, incidence),
          #* Discount rate assuming annual cycle length
          #* use 5% for LMIC from (Haacker, Hallett & Atun, 2020)
          discount_rates = c(0.05, 0.05, 0))

#* 3. Twig for 80% vaccination coverage
hpv_twig_80 <- twig() + 
  decisions(names = c(no_trt, FOV, GNV)) + 
  states(names = c(HV, HU, HPV_II, HPV_PI, PL, CL, Death),
         init_probs = c(0.8, 0.2, 0, 0, 0 , 0, 0),
         max_cycles = c(1, 1, 1, 1, 1, 1, 1)) + 
  event(name = death_event,
        options = c(yes, none),
        probs = c(p_die, leftover), 
        transitions = c(Death, second_event)) +
  event(name = second_event,
        options     = c(recover,   get_II,   get_PI,   get_PL,   get_CL,   none),
        probs       = c(p_recover, p_get_II, p_get_PI, p_get_PL, p_get_CL, leftover),
        transitions = c(HV,        HPV_II,   HPV_PI,   PL,       CL,       stay)) +
  payoffs(names = c(cost, utility, incidence),
          #* Discount rate assuming annual cycle length
          #* use 5% for LMIC from (Haacker, Hallett & Atun, 2020)
          discount_rates = c(0.05, 0.05, 0))


# 05 Run twigs ------------------------------------------------------------


results_40 <- run_twig(twig_obj = hpv_twig_40, 
                       params = df_params_40, 
                       n_cycles = 94, 
                       verbose = TRUE)

results_60 <- run_twig(twig_obj = hpv_twig_60, 
                       params = df_params_60, 
                       n_cycles = 94, 
                       verbose = TRUE)

results_80 <- run_twig(twig_obj = hpv_twig_80, 
                       params = df_params_80, 
                       n_cycles = 94, 
                       verbose = TRUE)



# 06 Generate ICERs -------------------------------------------------------

# Convert into data frame
df_payoffs_40 <- as.data.frame(results_40$sim_ev)
df_payoffs_60 <- as.data.frame(results_60$sim_ev)
df_payoffs_80 <- as.data.frame(results_80$sim_ev)

df_payoffs_40 <- tibble::rownames_to_column(df_payoffs_40, var = "strategy")
df_payoffs_60 <- tibble::rownames_to_column(df_payoffs_60, var = "strategy")
df_payoffs_80 <- tibble::rownames_to_column(df_payoffs_80, var = "strategy")

df_payoffs_all <- bind_rows(`40` = df_payoffs_40,
                            `60` = df_payoffs_60,
                            `80` = df_payoffs_80,
                            .id = "coverage")

df_payoffs_filt <- df_payoffs_all %>% 
  filter(strategy != "no_trt" | !(coverage %in% c("60", "80"))) %>% 
  mutate(str_name = stringr::str_c(strategy, coverage, sep = "_"))


df_payoffs_filt$str_name[1] <- "no_trt"
# df_payoffs_filt$str_name[6] <- "FOV_80" # Standard of Care

df_payoffs_filt <- rownames_to_column(df_payoffs_filt, var = "decision")

# df_payoffs_filt <- dplyr::rename(df_payoffs_filt, decision = strategy)


df_ICERs <- twig::calculate_icers(df_payoffs_filt)

#* Use the same names as the output of `dampack::calculate_icers()`
colnames(df_ICERs) <- c("Strategy",
                        "Cost",
                        "Effect",
                        "Inc_Cost",
                        "Inc_Effect",
                        "ICER",
                        "Status")

df_ICERs <- left_join(x = df_ICERs,
                      y = select(df_payoffs_filt, decision, str_name),
                      by = join_by(Strategy == decision)) %>% 
  select(-Strategy) %>% 
  rename(Strategy = str_name)


plt_ICERs <- darthtools::plot_icers(df_ICERs, label = "all") + 
  theme_bw(base_size = 18) + 
  theme(legend.position = "bottom")




# 08 Plot cohort trace ----------------------------------------------------

m_M_no_trt_40 <- results_40$markov_trace[,, "no_trt"]
m_M_FOV_40    <- results_40$markov_trace[,, "FOV"]
m_M_GNV_40    <- results_40$markov_trace[,, "GNV"]
m_M_no_trt_60 <- results_60$markov_trace[,, "no_trt"]
m_M_FOV_60    <- results_60$markov_trace[,, "FOV"]
m_M_GNV_60    <- results_60$markov_trace[,, "GNV"]
m_M_no_trt_80 <- results_80$markov_trace[,, "no_trt"]
m_M_FOV_80    <- results_80$markov_trace[,, "FOV"]
m_M_GNV_80    <- results_80$markov_trace[,, "GNV"]

# Merge HV and HU into HU when "no trt"
m_M_no_trt_40[, "HU"] <- rowSums(m_M_no_trt_40[, c("HV", "HU")])
m_M_no_trt_60[, "HU"] <- rowSums(m_M_no_trt_60[, c("HV", "HU")])
m_M_no_trt_80[, "HU"] <- rowSums(m_M_no_trt_80[, c("HV", "HU")])

# Remove HV and HU columns
m_M_no_trt_40[, "HV"] <- NA
m_M_no_trt_60[, "HV"] <- NA
m_M_no_trt_80[, "HV"] <- NA

# Define vector of colum names
v_names_no_trt <- c("HV", "HU", "HPV-II", "HPV-PI", "PL", "CL", "Death")
v_names_vacc   <- c("HV", "HU", "HPV-II", "HPV-PI", "PL", "CL", "Death")

# Update column names
colnames(m_M_no_trt_40) <- v_names_vacc
colnames(m_M_no_trt_60) <- v_names_vacc
colnames(m_M_no_trt_80) <- v_names_vacc
colnames(m_M_FOV_40)    <- v_names_vacc 
colnames(m_M_GNV_40)    <- v_names_vacc 
colnames(m_M_FOV_60)    <- v_names_vacc 
colnames(m_M_GNV_60)    <- v_names_vacc 
colnames(m_M_FOV_80)    <- v_names_vacc 
colnames(m_M_GNV_80)    <- v_names_vacc


df_M_all <- bind_rows(no_trt_40 = as.data.frame(m_M_no_trt_40),
                      FOV_40    = as.data.frame(m_M_FOV_40),
                      GNV_40    = as.data.frame(m_M_GNV_40),
                      no_trt_60 = as.data.frame(m_M_no_trt_60),
                      FOV_60    = as.data.frame(m_M_FOV_60),
                      GNV_60    = as.data.frame(m_M_GNV_60),
                      no_trt_80 = as.data.frame(m_M_no_trt_80),
                      FOV_80    = as.data.frame(m_M_FOV_80),
                      GNV_80    = as.data.frame(m_M_GNV_80),
                      .id = "source")

# Add format to dataframe
df_M_all <- cbind(cycle = rep(0:94, 9), df_M_all)
rownames(df_M_all) <- 1:nrow(df_M_all)





# Define vectors of colors and linetypes associated to specific health states
v_color_values <- c("HV"     = "#1B9E77",
                    "HU"     = "#D95F02", 
                    "HPV-II" = "#7570B3",
                    "HPV-PI" = "#E7298A", 
                    "PL"     = "#66A61E",
                    "CL"     = "#E6AB02",
                    "Death"  = "#A6761D")

v_linetype_values <- c("HV"     = "solid",
                       "HU"     = "22", 
                       "HPV-II" = "42",
                       "HPV-PI" = "44", 
                       "PL"     = "13",
                       "CL"     = "1343",
                       "Death"  = "73")


#* *40%* Vaccination coverage
plt_no_trt_40 <- plot_m_M(m_M_no_trt_40) +
  scale_color_manual(values = v_color_values) +
  scale_linetype_manual(values = v_linetype_values) + 
  theme(legend.position = "none") + 
  labs(subtitle = "No vaccination")

plt_FOV_40 <- plot_m_M(m_M_FOV_40) +
  scale_color_manual(values = v_color_values) +
  scale_linetype_manual(values = v_linetype_values) +
  theme(legend.position = "none") + 
  labs(subtitle = "Female-only vaccination")

plt_GNV_40 <- plot_m_M(m_M_GNV_40) +
  scale_color_manual(values = v_color_values) +
  scale_linetype_manual(values = v_linetype_values) +
  labs(subtitle = "Gender-neutral vaccination")

#* *60%* Vaccination coverage
plt_no_trt_60 <- plot_m_M(m_M_no_trt_60) +
  scale_color_manual(values = v_color_values) +
  scale_linetype_manual(values = v_linetype_values) + 
  theme(legend.position = "none") + 
  labs(subtitle = "No vaccination")

plt_FOV_60 <- plot_m_M(m_M_FOV_60) +
  scale_color_manual(values = v_color_values) +
  scale_linetype_manual(values = v_linetype_values) +
  theme(legend.position = "none") + 
  labs(subtitle = "Female-only vaccination")

plt_GNV_60 <- plot_m_M(m_M_GNV_60) +
  scale_color_manual(values = v_color_values) +
  scale_linetype_manual(values = v_linetype_values) +
  labs(subtitle = "Gender-neutral vaccination")

#* *80%* Vaccination coverage
plt_no_trt_80 <- plot_m_M(m_M_no_trt_80) +
  scale_color_manual(values = v_color_values) +
  scale_linetype_manual(values = v_linetype_values) + 
  theme(legend.position = "none") + 
  labs(subtitle = "No vaccination")

plt_FOV_80 <- plot_m_M(m_M_FOV_80) +
  scale_color_manual(values = v_color_values) +
  scale_linetype_manual(values = v_linetype_values) +
  theme(legend.position = "none") + 
  labs(subtitle = "Female-only vaccination")

plt_GNV_80 <- plot_m_M(m_M_GNV_80) +
  scale_color_manual(values = v_color_values) +
  scale_linetype_manual(values = v_linetype_values) +
  labs(subtitle = "Gender-neutral vaccination")


# Merge plots into a patchworks per % coverage
plt_M_40 <- plt_no_trt_40/plt_FOV_40/plt_GNV_40 +
  patchwork::plot_annotation(tag_levels = "A", tag_suffix = ")")

plt_M_60 <- plt_no_trt_60/plt_FOV_60/plt_GNV_60 +
  patchwork::plot_annotation(tag_levels = "A", tag_suffix = ")")

plt_M_80 <- plt_no_trt_80/plt_FOV_80/plt_GNV_80 +
  patchwork::plot_annotation(tag_levels = "A", tag_suffix = ")")





# 09 PSA ------------------------------------------------------------------


## 09.01 Define general parameters -----------------------------------------

# Define number of samples to be used in PSA
n_samp <- 5000


# Define model parameters using calibrated values
df_params_PSA <- dplyr::tibble(
  
  # # Proportion of vaccinated people by age
  # prop_vac_women = 0.8, #* *Assumption*
  # prop_vac_men   = 0.4, #* *Assumption*
  
  # Transition parameters parameters
  r_H_II     = 0.016260416,  #* *calibrated value* check script `analysis/02.01_calibration.R`
  r_II_PI    = 0.005949991,  #* *calibrated value* check script `analysis/02.01_calibration.R`
  r_PI_PL    = 0.010190967,  #* *calibrated value* check script `analysis/02.01_calibration.R`
  r_PL_CL    = 0.006369672,  #* *calibrated value* check script `analysis/02.01_calibration.R`
  cc_c_mort_rate = 6.58/1e5, # CC crude mortality rate. Source: section 3.3.3 of Human Papillomavirus and Related Diseases Report MEXICO
  hr_PL      = 1.1,  #* *assumed values* HR applied to baseline mort rate when in PL
  hr_CL      = 1.3,  #* *assumed values* HR applied to baseline mort rate when in CL
  hr_hpv_vac = 0.05, #* *PSA* HR representing protective effect of vaccination among females (assuming 95% of efficacy, Baussano, et al., 2017)
  hr_hpv_men = 0.2,  #* *PSA* HR representing protective HPV herd effect of gender neutral vaccination
  
  # Recovery rates
  r_II_HV = 0.7, # source: https://pmc.ncbi.nlm.nih.gov/articles/PMC6223532/#:~:text=The%20lrHPVs%20are%20mainly%20associated,20%25%20of%20cervical%20cancers%20worldwide.&text=Approximately%2070%25%20of%20HPV%20infections,in%20the%20remainder%20of%20cases.&text=Although%20persistent%20infection%20with%20hrHPV,effective%20cell%2Dmediated%20immune%20response.&text=Therefore%2C%20HIV%2Dpositive%20individuals%20who,risk%20factors%20among%20Nigerian%20women.
  r_PI_HV = 0.5, #* *check* source: https://link.springer.com/article/10.1186/s12905-023-02764-8
  r_PL_HV = 0.4, #* *assumed values* (ordinal recovery rate)
  r_CL_HV = 0.2, #* *assumed values* (ordinal recovery rate)
  
  # Cost parameters (USD)
  #* *PSA* Use gamma distribution for PSA
  c_vacc = rgamma_mean_sd(n_samples = n_samp, mean_x = 80, sd_x = 36, N = 1),

  #* *PSA* cytology test values from table 1 of (https://www.scielo.org.mx/scielo.php?script=sci_arttext&pid=S0036-36342014000500017)
  c_pap_smear = rnorm(n = n_samp, mean = 14.045, sd = 1.43),
  c_HV    = c_vacc + c_pap_smear,     #* Vaccination + annual pap smear (Beal et al, 2014). June 2013 price converted to Jan 2026 assuming general CPI inflation
  c_HU    = c_pap_smear,              #* Assuming annual pap smear (Beal et al, 2014). June 2013 price converted to Jan 2026 assuming general CPI inflation
  
  #* *PSA* Cost of genital wart treatment, 
  #* assuming a gamma distribution with approx 95% CI between 250 to 500 USD
  c_genital_wart = rgamma_mean_sd(n_samples = n_samp, mean_x = 375 ,sd_x = 70, N = 1),
  c_PI   = c_genital_wart,  
  c_II   = c_genital_wart,  
  
  #* *PSA* Cost of CIN treatment
  #* assuming a gamma distribution with approx 95% CI between 2000 and 3000 USD
  c_PL   = rgamma_mean_sd(n_samples = n_samp, mean_x = 2500,sd_x = 255.1,N = 1), #* Cost of CIN treatment (Insinga, et al, 2007) MXN price from Jan 2005 to Jan 2026, then applied official MXN-USD exchange rate of Jan 2026
  
  #* *PSA* Assuming a gamma distribution with approx 95% CI between 5,000 to 13,000 USD (Granados-García et al., 2019; Granados-García et al., 2014)
  c_CL   = rnorm(n = n_samp, mean = 9000, sd = 2040.92), #* Cost of CC treatment - Total (Granados-García et al., 2019). June 2017 price converted to Jan 2026 assuming general CPI inflation
  c_D    = 0,                  #* Assuming the public healthcare provider
  
    
  # Utility parameters
  
  #* *PSA* utility weight when having genital warts
  #* Assuming a lognormal distribution with approx 95% CI between 0.89 and 0.97
  u_g_warts = rlnorm_disutility(n_samples = n_samp,  util_mean = 0.95,  util_sd = 0.15,  N = 50),
  
  ## Health states
  u_HV = 1,               # Time-dependent utility weight for healthy state 
  u_HU = 1,               # Time-dependent utility weight for healthy state
  u_II = u_g_warts,       # Utility weight for genital wart (u_g_warts)
  u_PI = u_g_warts,       # Utility weight for genital wart (u_g_warts)
  
  #* *PSA* utility weight when having Precancerous lesion (CIN 2/3)
  #* Assuming a lognormal distribution with approx 95% CI between 0.82 and 0.90
  u_PL = rlnorm_disutility(n_samples = n_samp,  util_mean = 0.87,  util_sd = 0.15,  N = 50),
  
  #* *PSA* utility weight when having cervical cancer
  #* Assuming a lognormal distribution with approx 95% CI between 0.6 and 0.8
  u_CL = rlnorm_disutility(n_samples = n_samp,  util_mean = 0.71,  util_sd = 0.35,  N = 50),  
)



## 08.02 Scenario-specific parameters --------------------------------------

df_params_PSA_40 <- df_params_PSA_60 <- df_params_PSA_80 <- df_params_PSA

# assuming 40% coverage and 95% of efficacy
# Value calculated using values from Table 2 of (Baussano, et al., 2017)
df_params_PSA_40$hr_hpv_men <- 1/(0.67/(0.4*0.95)) 

# assuming 60% coverage and 95% of efficacy
# Value calculated using values from Table 2 of (Baussano, et al., 2017)
df_params_PSA_60$hr_hpv_men <- 1/(0.90/(0.6*0.95)) 

# assuming 80% coverage and 95% of efficacy
# Value calculated using values from Table 2 of (Baussano, et al., 2017)
df_params_PSA_80$hr_hpv_men <- 1/(1.00/(0.8*0.95)) 


## 08.03 Run PSA models ----------------------------------------------------


# 40% vaccination coverage
results_PSA_40 <- run_twig(twig_obj = hpv_twig_40, 
                           params = df_params_PSA_40, 
                           n_cycles = 94)

# 60% vaccination coverage
results_PSA_60 <- run_twig(twig_obj = hpv_twig_60, 
                           params = df_params_PSA_60, 
                           n_cycles = 94)

# 80% vaccination coverage
results_PSA_80 <- run_twig(twig_obj = hpv_twig_80, 
                           params = df_params_PSA_80, 
                           n_cycles = 94)



# 10 PSA plots ------------------------------------------------------------


# results_PSA_40 <- readRDS(file = "output/twig_results_PSA_40.rds")
# results_PSA_60 <- readRDS(file = "output/twig_results_PSA_60.rds")
# results_PSA_80 <- readRDS(file = "output/twig_results_PSA_80.rds")

## 10.01 CEAC --------------------------------------------------------------

# Cost effectiveness Acceptability curve

## 40% Vaccination coverage
plt_CEAC_40 <- plot_ceac(results_PSA_40$sim_ev, 
          wtp_range = seq(0, 100000, by = 1000)) + 
  theme(legend.position = "bottom")

## 60% Vaccination coverage
plt_CEAC_60 <- plot_ceac(results_PSA_60$sim_ev, 
          wtp_range = seq(0, 100000, by = 1000)) + 
  theme(legend.position = "bottom")

## 80% Vaccination coverage
plt_CEAC_80 <- plot_ceac(results_PSA_80$sim_ev, 
          wtp_range = seq(0, 100000, by = 1000)) + 
  theme(legend.position = "bottom")



## 10.02 CE scatterplot ----------------------------------------------------

plt_ce_scatter_40 <- plot_ce_scatterplot(l_out_PSA = results_PSA_40)
plt_ce_scatter_60 <- plot_ce_scatterplot(l_out_PSA = results_PSA_60)
plt_ce_scatter_80 <- plot_ce_scatterplot(l_out_PSA = results_PSA_80)

# 10 Save outputs ---------------------------------------------------------


## 10.01 Data --------------------------------------------------------------

## Twig outputs
saveRDS(object = results_40, file = "output/twig_results_40.rds")
saveRDS(object = results_60, file = "output/twig_results_60.rds")
saveRDS(object = results_80, file = "output/twig_results_80.rds")

## Cohort traces
saveRDS(object = df_M_all, file = "output/df_M_all.rds")

# ICERs
saveRDS(object = df_ICERs, file = "output/df_ICERs.rds")
write.csv(x = df_ICERs, file = "output/df_ICERs.csv")

# Save PSA results
saveRDS(object = results_PSA_40, file = "output/twig_results_PSA_40.rds")
saveRDS(object = results_PSA_60, file = "output/twig_results_PSA_60.rds")
saveRDS(object = results_PSA_80, file = "output/twig_results_PSA_80.rds")

## 10.02 Plots -------------------------------------------------------------

## Cost Effectiveness efficient frontier
ggsave(plot = plt_ICERs, filename = "figs/plt_ICERs.png", width = 10, height = 6)

## Cohort traces
ggsave(plot = plt_M_40, filename = "figs/plt_M_40.png", width = 8, height = 12)
ggsave(plot = plt_M_60, filename = "figs/plt_M_60.png", width = 8, height = 12)
ggsave(plot = plt_M_80, filename = "figs/plt_M_80.png", width = 8, height = 12)

## CEAC
ggsave(plot = plt_CEAC_40, filename = "figs/plt_CEAC_40.png", width = 10, height = 6)
ggsave(plot = plt_CEAC_60, filename = "figs/plt_CEAC_60.png", width = 10, height = 6)
ggsave(plot = plt_CEAC_80, filename = "figs/plt_CEAC_80.png", width = 10, height = 6)

## CE scatterplots
ggsave(plot = plt_ce_scatter_40, filename = "figs/plt_ce_scatter_40.png", width = 10, height = 6)
ggsave(plot = plt_ce_scatter_60, filename = "figs/plt_ce_scatter_60.png", width = 10, height = 6)
ggsave(plot = plt_ce_scatter_80, filename = "figs/plt_ce_scatter_80.png", width = 10, height = 6)



#* Notes for the future:
#* 
#* In the CE scatterplot
#* -  How many are dominated?
#* - How many are saving money?
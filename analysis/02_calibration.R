#* Title: 01_decision_twig
#* 
#* Code function: 
#*    This script contains the code used to calibrate the decision twig to 
#*    HPV epidemiological outcomes in Mexico.
#*    
#*    The females in the model will have a *minimal age of 16*, assuming they 
#*    could already have received the HPV vaccine.
#*    
#*    - Discount rate: 5% for countries in development (Gómez et al., 2016; Haacker, et al. 2020) 
#* 
#* Creation date: April 9th, 2026
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


#* The following functions are defined in `R/twig_funs.R`
#*    - defining the transition rates,
#*    - defining the costs,
#*    - defining the utilities,
#*    - calculating the incidence,
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


# 02.02 HPV prevalence in Mexico ------------------------------------------

#* *HPV prevalence in Mexico*
#* source: Appendix of https://www.sciencedirect.com/science/article/pii/S2667193X25001668?via%3Dihub
df_HPV_prev <- read_csv(file = "data-raw/HPV_overall_prev_mx_2023.csv")

# Get vector of HPV prevalence by age
v_HPV_prev <- rep(x = df_HPV_prev$prevalence,
                  times = c(9, 10, 10, 10, 10, 10, 35))

names(v_HPV_prev) <- 16:109

#* Calculate overall prevalence
#* - Calculation matches with statement of "Prevalence of circulating HPV genotypes" section of https://www.sciencedirect.com/science/article/pii/S2667193X25001668?via%3Dihub
df_HPV_overall <- df_HPV_prev %>% 
  mutate(n_cases = n*prevalence) %>% 
  summarise(`Age group` = "All",
            n_pop = sum(n),
            n_cases = sum(n_cases)) %>% 
  mutate(prevalence = n_cases/n_pop)

n_overall_HPV_prev <- df_HPV_overall$prevalence

## 02.03 HPV PI by age -----------------------------------------------------

#* *Prevalence of Persistent HPV Infection by age*
#* source: Table 3 of https://pmc.ncbi.nlm.nih.gov/articles/PMC11519437/#:~:text=The%20persistent%20HPV%20infection%20at,potentially%20allow%20HPV%20to%20evade.
df_HPV_PI_prev <- data.frame(age = c("≤30", "31-40", "41-50", ">50"),
                             prevalence = c(0.2393, 0.2122, 0.2113, 0.3255))

# Get vector of persistent HPV infection by age
v_HPV_PI_prev <- rep(x = df_HPV_PI_prev$prevalence,
                     times = c(15, 10, 10, 59))

names(v_HPV_PI_prev) <- 16:109



## 02.04 Utility weights in healthy state by age ---------------------------

#* Source: appendix table 7 of Elbasha et al, 2007
#* https://pubmed.ncbi.nlm.nih.gov/17370513/

v_u_H_0 <- c(0.93, 0.91, 0.89, 0.86, 0.80, 0.78, 0.70)

names(v_u_H_0) <- c("16-17", "18-34", "35-44", "45-54", "55-64", "65-74", "75+")

v_u_H <- rep(x = v_u_H_0, times = c(2, 17, 10, 10, 10, 10, 35))

names(v_u_H) <- 16:109


## 02.05 Population by sex in Mexico ---------------------------------------

# #* Population at the start of the year
# #* Source: 
# #* https://www.gob.mx/conapo/documentos/bases-de-datos-de-la-conciliacion-demografica-1950-a-2019-y-proyecciones-de-la-poblacion-de-mexico-2020-a-2070?idiom=es
# df_pob_raw  <- readxl::read_xlsx(path = "data-raw/ConDem50a19_ProyPob20a70/0_Pob_Inicio_1950_2070.xlsx")
# 
# #* Prepare and filter data
# #* - translate column names to english
# #* - keep population in calendar years from 2020 to 2030
# #* - keep only individuals aged 16
# #* - data for the whole country
# df_pop_2026 <- df_pob_raw %>% 
#   rename(row = RENGLON,
#          year = AÑO,
#          state = ENTIDAD,
#          geo_code = CVE_GEO,
#          age = EDAD,
#          sex = SEXO, 
#          population = POBLACION) %>% 
#   filter(year %in% 2020:2030,
#          state == "República Mexicana",
#          age == 16)
# 
# # Calculate the mean population by sex between 2020 and 2030
# df_pop_2026_sum <- df_pop_2026 %>% 
#   group_by(state, age, sex) %>% 
#   summarise(pop = mean(population))
# 
# # Extract number of men and women aged 16 in Mexico
# n_men <- df_pop_2026_sum %>% 
#   filter(sex == "Hombres") %>% 
#   pull(pop)
# 
# n_women <- df_pop_2026_sum %>% 
#   filter(sex == "Mujeres") %>% 
#   pull(pop)


## 02.06 CC mortality ------------------------------------------------------

#* Source of information: 
#* Section 3.3.3 of Human Papillomavirus and Related Diseases Report MEXICO
#* - https://hpvcentre.net/datastatistics.php
#* - https://hpvcentre.net/statistics/reports/MEX.pdf?t=1775754116625

# CC Crude mortality rate
cc_c_mort_rate <- 6.58/1e5

# CC age-standardized mortality rate
cc_as_mort_rate <- 5.74/1e5

## CC cumulative risk (mortality) at age 75
cc_cum_mort_risk <- 0.63


# 02.07 CC incidence ------------------------------------------------------

#* Source of data
#* Section 3.3.1 of Human Papillomavirus and Related Diseases Report MEXICO
#* - https://hpvcentre.net/datastatistics.php
#* - https://hpvcentre.net/statistics/reports/MEX.pdf?t=1775754116625

# Cervical Cancer crude incidence rate in Mexico
cc_c_inc_rate <- 14.3

# 03 General parameters ---------------------------------------------------

# Number of simulations
n_sim <- 1e3

# Random samples of model parameters
df_params <- dplyr::tibble(
  
  # Average number of men and women aged 16 in Mexico between 2020 and 2030
  n_men = 1141675,
  n_women = 1106284,
  
  # # Proportion of vaccinated people by age
  # prop_vac_women = 0.8, #* *Assumption*
  # prop_vac_men   = 0.4, #* *Assumption*
  
  # Transition parameters parameters
  r_H_II     = 0.80, #* *made-up values*
  r_II_PI    = 0.50, #* *made-up values*
  r_PI_PL    = 0.4,  #* *made-up values*
  r_PL_CL    = 0.3,  #* *made-up values*
  cc_c_mort_rate = 6.58/1e5, # CC crude mortality rate. Source: section 3.3.3 of Human Papillomavirus and Related Diseases Report MEXICO
  hr_PL      = 1.1,  #* *made-up values* HR applied to baseline mort rate when in PL
  hr_CL      = 1.3,  #* *made-up values* HR applied to baseline mort rate when in CL
  hr_hpv_vac = 0.6,  #* HR representing protective effect of vaccination among females
  hr_hpv_men = 0.2,  #* HR representing protective HPV herd effect of vaccinating males
  
  # Recovery rates
  r_II_HV = 0.7, # source: https://pmc.ncbi.nlm.nih.gov/articles/PMC6223532/#:~:text=The%20lrHPVs%20are%20mainly%20associated,20%25%20of%20cervical%20cancers%20worldwide.&text=Approximately%2070%25%20of%20HPV%20infections,in%20the%20remainder%20of%20cases.&text=Although%20persistent%20infection%20with%20hrHPV,effective%20cell%2Dmediated%20immune%20response.&text=Therefore%2C%20HIV%2Dpositive%20individuals%20who,risk%20factors%20among%20Nigerian%20women.
  r_PI_HV = 0.5, #* *check* source: https://link.springer.com/article/10.1186/s12905-023-02764-8
  r_PL_HV = 0.4, #* *Made up* (ordinal recovery rate)
  r_CL_HV = 0.2, #* *Made up* (ordinal recovery rate)
  
  # Cost parameters (USD)
  #* https://www.saludpublica.mx/index.php/spm/article/view/7374
  #* https://pmc.ncbi.nlm.nih.gov/articles/PMC3025113/#Sec7
  #* https://www.scielosp.org/pdf/spm/2014.v56n5/502-510#:~:text=Objective.,of%20undetected%20cervical%20cancer%20cases.&text=cervical%20cancer;%20Mexico-,Granados%2DGarc%C3%ADa%20V%2C%20Flores%20YN%2C%20P%C3%A9rez%20R%2C%20Rudolph,;56:502%2D510.
  c_vacc = ((331*2)*1.13)/17.95, # Using cost of CENETEC 2023, considering inflation to Jan 2026 and MXN-USD exchange rate of Jan 2026
  # c_trt_FOV  = n_women*prop_vac_women*c_vacc,
  # c_trt_vacc = c_vacc*(n_women*prop_vac_women + n_men*prop_vac_men),
  
  c_HV    = c_vacc + 19.5,     #* Vaccination + annual pap smear (Beal et al, 2014). June 2013 price converted to Jan 2026 assuming general CPI inflation
  c_HU    = 19.5,              #* Assuming annual pap smear (Beal et al, 2014). June 2013 price converted to Jan 2026 assuming general CPI inflation
  c_PI   = (2000*2.52)/17.95,  #* Cost of genital wart treatment (Insinga, et al, 2007) MXN price from Jan 2005 to Jan 2026, then applied official MXN-USD exchange rate of Jan 2026
  c_II   = (2000*2.52)/17.95,  #* Cost of genital wart treatment (Insinga, et al, 2007) MXN price from Jan 2005 to Jan 2026, then applied official MXN-USD exchange rate of Jan 2026
  c_PL   = (15849*2.52)/17.95, #* Cost of CIN treatment (Insinga, et al, 2007) MXN price from Jan 2005 to Jan 2026, then applied official MXN-USD exchange rate of Jan 2026
  c_CL   = 952.54,             #* Mean monthly cost of CC treatment - all stages (Granados-García et al., 2019). June 2017 price converted to Jan 2026 assuming general CPI inflation
  c_D    = 0,                  #* Assuming the public healthcare provider
  
  # # Utility parameters
  u_HV = 1,               # Time-dependent utility weight for healthy state 
  u_HU = 1,               # Time-dependent utility weight for healthy state
  u_II = 0.91,            # Utility weight for genital wart (Elbasha et al. 2007)
  u_PI = 0.91,            # Utility weight for CIN 1 (Elbasha et al. 2007)
  u_PL = 0.87,            # Utility weight for CIN 2/3 (Elbasha et al. 2007)
  u_CL = (0.76 + 0.67)/2  # Average Utility weight for localized and regional cervical cancer (Elbasha et al. 2007)
)




# 04 Define twig ----------------------------------------------------------


hpv_twig <- twig() + 
  decisions(names = c(no_trt, FOV, GNV)) + 
  states(names = c(HV, HU, HPV_II, HPV_PI, PL, CL, Death),
         init_probs = c(0.8, 0.2, 0, 0, 0 , 0, 0),
         max_cycles = c(1, 1, 1, 1, 1, 1, 1)) + # cycle_in_state references the time spent in the tunnel
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



# 05 Calibration ----------------------------------------------------------

model_function <- function(p) {
  
  # Random samples of model parameters
  df_params <- dplyr::tibble(
    n_men = 1141675,   # Average number of men and women aged 16 in Mexico between 2020 and 2030
    n_women = 1106284, # Average number of men and women aged 16 in Mexico between 2020 and 2030
    
    # # Proportion of vaccinated people by age
    # prop_vac_women = 0.8, #* *Assumption*
    # prop_vac_men   = 0.4, #* *Assumption*
    
    # Transition parameters parameters
    r_H_II         = p[1],  #* *calibrated values*
    r_II_PI        = p[2],  #* *calibrated values*
    r_PI_PL        = p[3],  #* *calibrated values*
    r_PL_CL        = p[4],  #* *calibrated values*
    cc_c_mort_rate = 6.58/1e5, # CC crude mortality rate. Source: section 3.3.3 of Human Papillomavirus and Related Diseases Report MEXICO
    hr_PL      = 1.1,  #* *made-up values* HR applied to baseline mort rate when in PL
    hr_CL      = 1.3,  #* *made-up values* HR applied to baseline mort rate when in CL
    hr_hpv_vac = 0.2,  #* HR representing protective effect of vaccination among females
    hr_hpv_men = 0.6,  #* HR representing protective HPV herd effect of vaccinating males
    
    # Recovery rates
    r_II_HV = 0.7, # source: https://pmc.ncbi.nlm.nih.gov/articles/PMC6223532/#:~:text=The%20lrHPVs%20are%20mainly%20associated,20%25%20of%20cervical%20cancers%20worldwide.&text=Approximately%2070%25%20of%20HPV%20infections,in%20the%20remainder%20of%20cases.&text=Although%20persistent%20infection%20with%20hrHPV,effective%20cell%2Dmediated%20immune%20response.&text=Therefore%2C%20HIV%2Dpositive%20individuals%20who,risk%20factors%20among%20Nigerian%20women.
    r_PI_HV = 0.5, #* *check* source: https://link.springer.com/article/10.1186/s12905-023-02764-8
    r_PL_HV = 0.4, #* *Made up* (ordinal recovery rate)
    r_CL_HV = 0.2, #* *Made up* (ordinal recovery rate)
    
    # Cost parameters (USD)
    #* https://www.saludpublica.mx/index.php/spm/article/view/7374
    #* https://pmc.ncbi.nlm.nih.gov/articles/PMC3025113/#Sec7
    #* https://www.scielosp.org/pdf/spm/2014.v56n5/502-510#:~:text=Objective.,of%20undetected%20cervical%20cancer%20cases.&text=cervical%20cancer;%20Mexico-,Granados%2DGarc%C3%ADa%20V%2C%20Flores%20YN%2C%20P%C3%A9rez%20R%2C%20Rudolph,;56:502%2D510.
    c_vacc = ((331*2)*1.13)/17.95, # Using cost of CENETEC 2023, considering inflation to Jan 2026 and MXN-USD exchange rate of Jan 2026
    # c_trt_FOV  = n_women*prop_vac_women*c_vacc,
    # c_trt_vacc = c_vacc*(n_women*prop_vac_women + n_men*prop_vac_men),
    
    c_HV    = 41.67465 + 19.5,     #* Vaccination + annual pap smear (Beal et al, 2014). June 2013 price converted to Jan 2026 assuming general CPI inflation
    c_HU    = 19.5,              #* Assuming annual pap smear (Beal et al, 2014). June 2013 price converted to Jan 2026 assuming general CPI inflation
    c_PI   = (2000*2.52)/17.95,  #* Cost of genital wart treatment (Insinga, et al, 2007) MXN price from Jan 2005 to Jan 2026, then applied official MXN-USD exchange rate of Jan 2026
    c_II   = (2000*2.52)/17.95,  #* Cost of genital wart treatment (Insinga, et al, 2007) MXN price from Jan 2005 to Jan 2026, then applied official MXN-USD exchange rate of Jan 2026
    c_PL   = (15849*2.52)/17.95, #* Cost of CIN treatment (Insinga, et al, 2007) MXN price from Jan 2005 to Jan 2026, then applied official MXN-USD exchange rate of Jan 2026
    c_CL   = 952.54,             #* Mean monthly cost of CC treatment - all stages (Granados-García et al., 2019). June 2017 price converted to Jan 2026 assuming general CPI inflation
    c_D    = 0,                  #* Assuming the public healthcare provider
    
    # # Utility parameters
    u_HV = 1,               # Time-dependent utility weight for healthy state 
    u_HU = 1,               # Time-dependent utility weight for healthy state
    u_II = 0.91,            # Utility weight for genital wart (Elbasha et al. 2007)
    u_PI = 0.91,            # Utility weight for CIN 1 (Elbasha et al. 2007)
    u_PL = 0.87,            # Utility weight for CIN 2/3 (Elbasha et al. 2007)
    u_CL = (0.76 + 0.67)/2  # Average Utility weight for localized and regional cervical cancer (Elbasha et al. 2007)
  )
  
  results <- suppressMessages(suppressWarnings(
    run_twig(twig_obj = hpv_twig, 
             params = df_params, 
             n_cycles = 94, 
             verbose = TRUE)
  ))
  
  # Extract outputs to be calibrated
  
  ## Age-standardized HPV overall prevalence
  
  v_alive   <- (1 - results$markov_trace[,"Death", "FOV"])
  v_HPV_raw <- rowSums(results$markov_trace[,c("HPV_II", "HPV_PI"), "FOV"])
  n_HPV_prev_twig <- sum(v_HPV_raw)/sum(v_alive) # Proper 
  
  ## CC crude incidence rate (per 100,000)
  n_cc_inc_rate_twig <- results$sim_ev["FOV", "incidence"]*1e5
  
  l_out <- list(results      = results,
                HPV_prev     = n_HPV_prev_twig,
                CC_crude_inc = n_cc_inc_rate_twig)
  
  return(l_out)
}



gof <- function(p) {
  
  twig_out <- model_function(p)
  
  # compute goodness of fit as sum of squared errors
  # n_gof <- sum((pred_prob_age_25 - obs_prob_age_25)^2) 
  n_gof <- sum((twig_out$HPV_prev - n_overall_HPV_prev)^2,
               (twig_out$CC_crude_inc - cc_c_inc_rate)^2) 
  
  #return goodness of fit
  return(n_gof) 
}


# Check that function works
v_params <- c(r_H_II  = 0.80,
              r_II_PI = 0.50,
              r_PI_PL = 0.4,
              r_PL_CL = 0.3)

m_out <- model_function(p = v_params)

gof(p = v_params)

results <- m_out$results


# use optim function which returns prob of smoking and quitting that minimize the gof
optim_res <- optim(par = c(r_H_II  = 0.01,
                           r_II_PI = 0.01,
                           r_PI_PL = 0.01,
                           r_PL_CL = 0.01), 
                   fn = gof)

# optim_res$par
# r_H_II     r_II_PI     r_PI_PL     r_PL_CL 
# 0.016260416 0.005949991 0.010190967 0.006369672 

v_params_calib <- c(r_H_II = 0.016260416,
                    r_II_PI = 0.005949991,
                    r_PI_PL = 0.010190967,
                    r_PL_CL = 0.006369672)

# calibrated_results <- model_function(p = optim_res$par)
calibrated_results <- model_function(p = v_params_calib)


# 06 Save outputs ---------------------------------------------------------

saveRDS(object = optim_res, 
        file = "output/optim_res.rds")

saveRDS(object = calibrated_results, 
        file = "output/l_calibrated_results.rds")

saveRDS(object = optim_res$par,
        file = "output/l_calibrated_params.rds")



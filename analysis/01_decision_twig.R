#* Title: 01_decision_twig
#* 
#* Code function: 
#*    This script contains the code of the initial version of the decision
#*    twig used to model the effect of male HPV vaccination in Mexico
#*    
#*    The females in the model will have a *minimal age of 16*, assuming they 
#*    could already have received the HPV vaccine.
#*    
#*    - Discount rate: 5% for countries in development (Gómez et al., 2016; Haacker, et al. 2020) 
#* 
#* Creation date: March 16 2026
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

n_HPV_prev <- df_HPV_overall$prevalence

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

#* Population at the start of the year
#* Source: 
#* https://www.gob.mx/conapo/documentos/bases-de-datos-de-la-conciliacion-demografica-1950-a-2019-y-proyecciones-de-la-poblacion-de-mexico-2020-a-2070?idiom=es
df_pob_raw  <- readxl::read_xlsx(path = "data-raw/ConDem50a19_ProyPob20a70/0_Pob_Inicio_1950_2070.xlsx")

#* Prepare and filter data
#* - translate column names to english
#* - keep population in calendar years from 2020 to 2030
#* - keep only individuals aged 16
#* - data for the whole country
df_pop_2026 <- df_pob_raw %>% 
  rename(row = RENGLON,
         year = AÑO,
         state = ENTIDAD,
         geo_code = CVE_GEO,
         age = EDAD,
         sex = SEXO, 
         population = POBLACION) %>% 
  filter(year %in% 2020:2030,
         state == "República Mexicana",
         age == 16)

# Calculate the mean population by sex between 2020 and 2030
df_pop_2026_sum <- df_pop_2026 %>% 
  group_by(state, age, sex) %>% 
  summarise(pop = mean(population))

# Extract number of men and women aged 16 in Mexico
n_men <- df_pop_2026_sum %>% 
  filter(sex == "Hombres") %>% 
  pull(pop)

n_women <- df_pop_2026_sum %>% 
  filter(sex == "Mujeres") %>% 
  pull(pop)

# 03 General parameters ---------------------------------------------------



# Number of simulations
n_sim <- 1e3

# Random samples of model parameters
df_params <- dplyr::tibble(
  
  # Average number of men and women aged 16 in Mexico between 2020 and 2030
  n_men = n_men,
  n_women = n_women,
  
  # # Proportion of vaccinated people by age
  # prop_vac_women = 0.8, #* *Assumption*
  # prop_vac_men   = 0.4, #* *Assumption*
  
  # Transition parameters parameters
  r_H_II     = 0.22, #* *made-up values*
  r_II_PI    = 0.26, #* *made-up values*
  r_PI_PL    = 0.3,  #* *made-up values*
  r_PL_CL    = 0.4,  #* *made-up values*
  hr_CL      = 1.1,  #* HR applied to baseline mort rate when in CL
  hr_PL      = 1.1,  #* HR applied to baseline mort rate when in PL
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
  # c_trt_SoQ  = n_women*prop_vac_women*c_vacc,
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
  decisions(names = c(SoQ, male_vaccination)) + 
  states(names = c(HV, HU, HPV_II, HPV_PI, PL, CL, Death),
         init_probs = c(0.8, 0.2, 0, 0, 0 , 0, 0),
         max_cycles = c(1, 1, 1, 1, 1, 1, 1)) + # cycle_in_state references the time spent in the tunnel
  event(name = death_event,
        options = c(yes, none),
        probs = c(p_die, leftover), 
        transitions = c(Death, second_event)) +
  event(name = second_event,
        options     = c(recover, get_II, get_PI, get_PL, get_CL, none),
        probs       = c(p_recover, p_get_II, p_get_PI, p_get_PL, p_get_CL, leftover),
        transitions = c(HV, HPV_II, HPV_PI, PL, CL, stay)) +
  payoffs(names = c(cost, utility),
          #* Discount rate assuming annual cycle length
          #* use 5% for LMIC from (Haacker, Hallett & Atun, 2020)
          discount_rates = c(0.05, 0.05)) 


# 05 Define twig functions ------------------------------------------------

## 05.01 Transition functions ----------------------------------------------

p_die <- function(state, cycle, hr_PL, hr_CL) {
  
  r_mort_0 <- v_r_mort[cycle]
  
  # Calculate case-specific mortality rates
  r_mort_1 <- r_mort_0 *  # Baseline mortality if no (pre)cancerous lesion
    hr_PL^(state == "PL") * # Mortality rate when in PL
    hr_CL^(state == "CL") # Mortality rate when in CL
  
  # Convert rates into probabilities
  rate2prob(r_mort_1)
}


#* Function defining the probability of recovery
#* 
#* *Assumptions*: 
#* 1. current: The recovery rate is constant over time
#*     - (consider applying) In PI there is a 15% decrease in the clearance rate for every 5-year increase in age  https://www.sciencedirect.com/science/article/pii/S0090825819312703
#* 2. There is an *ordinal recovery rate* (i.e., higher stages -> lower recovery rates)

p_recover <- function(state, r_II_HV, r_PI_HV, r_PL_HV, r_CL_HV) {
  
  
  r_recovery <- r_II_HV*(state == "II") +
    r_PI_HV*(state == "PI") +
    r_PL_HV*(state == "PL") +
    r_CL_HV*(state == "CL")
  
  # Convert rates into probabilities
  rate2prob(r_recovery)
}


#* Function defining the probability of obtaining Incident HPV infection (II)
p_get_II  <- function(state, cycle, decision, r_H_II, hr_hpv_vac, hr_hpv_men) {
  
  # Change protective effect based on decisions
  hr_hpv_men <- 1*hr_hpv_men^(decision == "male_vaccination")
  
  #* Define risk of geting HPV Incident Infection (II)
  r_H_II_1 <- r_H_II*hr_hpv_men*hr_hpv_vac*(state == "HV") +  # when vaccinated
    r_H_II*hr_hpv_men*(state == "HU")                         # when non vaccinated
  
  rate2prob(r_H_II_1)
}

#* Function defining the probability of obtaining Incident HPV infection (PI)
p_get_PI  <- function(state, cycle, r_II_PI) {
    rate2prob(r_II_PI)*(state == "II")
}

p_get_PL  <- function(state, cycle, r_PI_PL) {
  rate2prob(r_PI_PL)*(state == "PI")
}


#* Function defining the incidence rate of CC
p_get_CL  <- function(state, cycle, r_PL_CL) {
  rate2prob(r_PL_CL)*(state == "PL")
}




## 05.02 Cost and utility functions ----------------------------------------


# #* Empty cost and utility functions.
# #* - The model is running without being calibrated (March 30th, 2026)
# cost <- function(decision) {
#   return(c(0, 0))
# }
# 
# utility <- function(decision) {
#   return(c(0, 0))
# }


cost <- function(state, c_HV, c_HU, c_II, c_PI, c_PL, c_CL) {
  
    # #* Decision costs
  # #* Extra cost in vaccination (assuming a 50% sex proportion in the population)
  # c_decision <- c_trt_SoQ*(decision == "SoQ") + 
  #   c_trt_vacc*(decision == "male_vaccination")
  
  # State costs
  c_state <- c_HV*(state == "HV") +
    c_HU*(state == "HU") +
    c_II*(state == "II") +
    c_PI*(state == "PI") +
    c_PL*(state == "PL") +
    c_CL*(state == "CL")
    
  return(c_state)
}

utility <- function(state, u_HV, u_HU, u_II, u_PI, u_PL, u_CL) {
  
  # State utilities
  u_state <- u_HV*(state == "HV") +
    u_HU*(state == "HU") +
    u_II*(state == "II") +
    u_PI*(state == "PI") +
    u_PL*(state == "PL") +
    u_CL*(state == "CL")
  

  return(u_state)
}

results <- run_twig(twig_obj = hpv_twig, params = df_params, n_cycles = 20, verbose = TRUE)

#* *Questions*
#* - How can we use Twig to calibrate the model? (Check the code from class Nelder-Mead example)
#* 
#* How to incorporate decision costs?
#* - The cost of vaccinating the male population
#* 
#* How to incorporatetTime-dependent utility weight for healthy state
#* 
#* The most important stages in litterature are related to CC, CIN *update the model*?

#* *April 2nd Notes*
#* 1 _DONE_ Take average of regional-local CC to generate the costs of CL
#* 2 _DONE_ Take the average cost of treatments for genital warts
#* 4 _DONE_ Make the cost of the intervention at an individual level, then it can be multiplied by the population
#* 3 _DONE_ Update cost of the vaccine: 1 adjust for mexican inflation, 2 transform from MXN to USD
#* - Use the overall HPV prevalence as a calibration target for II, Use CC prevalence as a calibration target for CL
#* - First calibrate the model, then run a sensitivity analysis
#* - _DONE_ Assume a utility weight utility of 1 for Healthy state

#* Remember to use the published effects when vaccinating males and females

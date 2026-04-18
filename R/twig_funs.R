#* Title: twig_funs.R
#* 
#* Code function:
#*    The functions 
#*    - defining the transition rates,
#*    - defining the costs,
#*    - defining the utilities,
#*    - calculating the incidence,
#*    - plotting the cohort trace
#*    are defined in this script and will be loaded in multiple R scripts inside
#*    the `analysis` folder
#* 
#* Creation date: April 16 2026
#* Author: David U. Garibay Treviño


p_die <- function(state, cycle, hr_PL, hr_CL) {
  
  r_mort_0 <- v_r_mort[cycle]
  
  # Calculate case-specific mortality rates
  r_mort_1 <- r_mort_0*    # Baseline mortality if no (pre)cancerous lesion
    hr_PL^(state == "PL")* # Mortality rate when in PL
    hr_CL^(state == "CL")  # Mortality rate when in CL
  
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
  
  # Protective effect on women due to vaccination
  hr_hpv_vac <- hr_hpv_vac^(decision %in% c("GNV", "FOV"))
  # Protective herd effect on women due to Gender neutral vaccination
  hr_hpv_men <- hr_hpv_men^(decision == "GNV")
  
  #* Define risk of geting HPV Incident Infection (II)
  r_H_II_1 <- r_H_II*hr_hpv_men*hr_hpv_vac*(state == "HV") +  # when vaccinated
    r_H_II*hr_hpv_men*(state == "HU")                         # when non vaccinated
  
  rate2prob(r_H_II_1)
}

#* Function defining the probability of obtaining Incident HPV infection (PI)
p_get_PI  <- function(state, cycle, r_II_PI) {
  
  rate2prob(r_II_PI)*(state == "HPV_II")
}

p_get_PL  <- function(state, cycle, r_PI_PL) {
  rate2prob(r_PI_PL)*(state == "HPV_PI")
}


#* Function defining the incidence rate of CC
p_get_CL  <- function(state, cycle, r_PL_CL) {
  rate2prob(r_PL_CL)*(state == "PL")
}

# cost <- function(state, cycle, decision, c_HV, c_HU, c_II, c_PI, c_PL, c_CL, c_vacc) {
# 
#   # #* Decision costs
#   # #* Extra cost in vaccination (assuming a 50% sex proportion in the population)
#   # c_decision <- c_trt_FOV*(decision == "FOV") +
#   #   c_trt_vacc*(decision == "GNV")
# 
#   # Assume double cost of vaccination when applying Gender neutral vaccination
#   c_HV <- c_HU + c_vacc*(decision != "no_trt")*2^(decision == "GNV")
# 
#   # State costs
#   c_state <- c_HV*(state == "HV") +
#     c_HU*(state == "HU") +
#     c_II*(state == "HPV_II") +
#     c_PI*(state == "HPV_PI") +
#     c_PL*(state == "PL") +
#     c_CL*(state == "CL")
# 
#   return(c_state)
# }

cost <- function(state, cycle, decision, c_HV, c_HU, c_II, c_PI, c_PL, c_CL, c_vacc) {

  # Decision costs

  #* How the decision cost (vaccination) is applied:
  #* - The cost is only applied in cycle 1,
  #* - The cost only applies to those who are in the HV state,
  #* - The cost does not apply if the decision is no treatment at all
  #* - If there is GNV, then the cost is doubled (assuming a 50% proportion by sex in the pop)
  c_decision <- (cycle == 1)*(state == "HV")*c_vacc*(decision != "no_trt")*2^(decision == "GNV")

  # State costs
  c_state <- c_HV*(state == "HV") +
    c_HU*(state == "HU") +
    c_II*(state == "HPV_II") +
    c_PI*(state == "HPV_PI") +
    c_PL*(state == "PL") +
    c_CL*(state == "CL")

  return(c_decision + c_state)
}


utility <- function(state, u_HV, u_HU, u_II, u_PI, u_PL, u_CL) {
  
  # State utilities
  u_state <- u_HV*(state == "HV") +
    u_HU*(state == "HU") +
    u_II*(state == "HPV_II") +
    u_PI*(state == "HPV_PI") +
    u_PL*(state == "PL") +
    u_CL*(state == "CL")
  
  
  return(u_state)
}

incidence <- function(second_event) {
  
  1*(second_event == "get_CL")
  
  # Multiply per 100,000
}


# random samples from log-normal distribution
rlnorm_disutility <- function(n_samples, util_mean, util_sd, N){
  disutil_mean <- 1 - util_mean
  util_se <- util_sd / sqrt(N)
  adjustment_factor <- log((util_se/disutil_mean)^2+1)
  disutil_mean_adj <- log(disutil_mean) - 0.5 * adjustment_factor
  disutil_se_adj <- sqrt(adjustment_factor)
  1 - rlnorm(n = n_samples, meanlog = disutil_mean_adj, sdlog = disutil_se_adj)
}

# random samples from Gamma distribution 
rgamma_mean_sd <- function(n_samples, mean_x, sd_x, N){
  se_x <- sd_x / sqrt(N)
  shape <- (mean_x / se_x)^2
  scale <- se_x^2 / mean_x
  rgamma(n = n_samples, shape = shape, scale = scale)
}

# random samples from beta distribution 
rbeta_mean_sd <- function(n_samples, mean_x, N){
  rbeta(n = n_samples, shape1 = mean_x * N, shape2 = N * (1 - mean_x))
}

# random samples from lognormal risks
rlnorm_ci <- function(n_samples, mean_x, lb_x, ub_x){
  rlnorm(n_samples, meanlog = log(mean_x), 
         sdlog = (log(ub_x) - log(lb_x)) / (1.96 * 2))
}


#* Define personalized function to plot cohort trace
#* - this is a modified version of `darthtools::plot_trace()`
plot_m_M <- function(m_M) {
  
  df_M <- data.frame(Age = 16 + (0:(nrow(m_M) - 1)), 
                     m_M, 
                     check.names = F)
  
  v_names_states <- colnames(m_M)
  
  df_M_long <- tidyr::pivot_longer(df_M, 
                                   cols = -Age, 
                                   names_to = "Health State", 
                                   values_to = "value")
  
  df_M_long$`Health State` <- factor(df_M_long$`Health State`, 
                                     levels = v_names_states)
  
  gg_trace <- ggplot2::ggplot(df_M_long, aes(x = Age, 
                                             y = value, 
                                             color = `Health State`, 
                                             linetype = `Health State`)) + 
    ggplot2::geom_line(linewidth = 1) + 
    ggplot2::xlab("Age") + 
    ggplot2::ylab("Proportion") + 
    ggplot2::scale_x_continuous(breaks = darthtools::number_ticks(8)) + 
    ggplot2::scale_y_continuous(breaks = darthtools::number_ticks(5)) + 
    ggplot2::theme_bw(base_size = 14) + 
    ggplot2::theme(legend.position = "bottom", 
                   legend.background = element_rect(fill = NA)) +
    ggplot2::guides(color = guide_legend(nrow = 1))
  
  return(gg_trace)
}


plot_ce_scatterplot <- function(l_out_PSA) {
  
  # Scatterplot
  df <- as.data.frame(as.table(l_out_PSA$sim_ev)) %>%
    setNames(c("decision", "payoff", "sim", "value")) %>%
    pivot_wider(names_from = payoff, values_from = value)
  
  head(df)
  # Compute incremental outcomes vs. NoTreatment
  df_inc <- df %>%
    group_by(sim) %>%
    mutate(
      cost_NoTx    = cost[decision == "no_trt"],
      util_NoTx    = utility[decision == "no_trt"],
      dCost        = cost - cost_NoTx,
      dUtility     = utility - util_NoTx
    ) %>%
    ungroup() %>%
    filter(decision != "no_trt") %>%
    mutate(
      comparison = case_when(
        decision == "FOV" ~ "FOV vs No Treatment",
        decision == "GNV" ~ "GNV vs No Treatment"
      )
    )
  
  plt <- ggplot(df_inc, aes(x = dUtility, y = dCost, color = comparison)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      x = "Incremental Utility",
      y = "Incremental Cost",
      color = "Comparison"
      # title = "Cost-Effectiveness Plane"
    ) +
    theme_bw() + 
    theme(legend.position = "bottom")
  
  return(plt)
}

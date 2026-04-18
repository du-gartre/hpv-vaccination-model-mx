#* Title: 00_calc_annual_mort_prob
#* 
#* Code function:
#*    This script loads population proyections from Mexico. Specifically,
#*    it loads the population at the beginning of the year and the 
#*    number of deaths by year. The source of the original data is:
#*    
#*    https://www.gob.mx/conapo/documentos/bases-de-datos-de-la-conciliacion-demografica-1950-a-2019-y-proyecciones-de-la-poblacion-de-mexico-2020-a-2070?idiom=es
#*    
#*    With that information this script computes
#*    the cumulative annual mortality rate, the cumulative probability of death,
#*    and the annual probability of death.
#* 
#* Creation date: February 16 2026
#* Author: David U. Garibay Treviño


# 01 Initial Setup --------------------------------------------------------

## 01.01 Clean environment ------------------------------------------------
remove(list = ls())

#* Refresh environment memory
gc()

## 01.02 Load libraries ----------------------------------------------------
library(readr)
library(dplyr)
library(readxl)




# 02 Read data ------------------------------------------------------------

# Population at the start of the year
df_pob_raw  <- read_xlsx(path = "data-raw/ConDem50a19_ProyPob20a70/0_Pob_Inicio_1950_2070.xlsx")

# Registered deaths
df_mort_raw <- read_xlsx(path = "data-raw/ConDem50a19_ProyPob20a70/1_Defunciones_1950_2070.xlsx")



# 03 Merge data -----------------------------------------------------------


df_pob_mort_raw <- left_join(x = df_pob_raw,
                             y = df_mort_raw,
                             by = c("AÑO", "ENTIDAD", "CVE_GEO", "EDAD", "SEXO"))

# 04 Clean data -----------------------------------------------------------


df_pob_mort_1 <- df_pob_mort_raw %>% 
  # Remove RENGLON columns
  dplyr::select(- RENGLON.x, -RENGLON.y) %>% 
  # Only keep national info for females (mujeres) in year 2020
  dplyr::filter(ENTIDAD == "República Mexicana",
                SEXO == "Mujeres",
                AÑO == 2020) %>% 
  # Modify DEFUNCIONES to make sure that # deaths <= # population
  mutate(DEFUNCIONES = if_else(condition = DEFUNCIONES <= POBLACION,
                               true = DEFUNCIONES,
                               false = POBLACION)) %>% 
  # Calculate mortality rate
  mutate(mort_rate = DEFUNCIONES/POBLACION) %>% 
  # Rename column
  select(year = AÑO,
         state = ENTIDAD,
         age = EDAD,
         sex = SEXO,
         pop = POBLACION,
         deaths = DEFUNCIONES,
         mort_rate)
  
#* Calculate
#* - H_cum:  Cumulative annual mortality rate
#* - P_cum:  Cumulative annual probability of death
#* - p_inst: Annual probability of death
df_pob_mort_2 <- df_pob_mort_1 %>% 
  mutate(H_cum  = cumsum(mort_rate),
         P_cum  = 1 - exp(-H_cum),
         p_inst = c(P_cum[1], diff(P_cum)))


# Make sure everyone dies at age 109
df_pob_mort_2[df_pob_mort_2$age == 109,]$P_cum <- 1
df_pob_mort_2[df_pob_mort_2$age == 109,]$p_inst <- 1 - df_pob_mort_2[df_pob_mort_2$age == 108,]$P_cum

#* The sum of `p_inst` should be equal to 1.
#* If that does not happen this line will return an error
stopifnot(sum(df_pob_mort_2$p_inst) == 1)

# 04 Generate visualizations ----------------------------------------------

# Plot annual probability of death
plt_mort_prob <- ggplot(data = df_pob_mort_2,
                        mapping = aes(x = age, 
                                      y = p_inst)) +
  theme_bw(base_size = 18) + 
  geom_line(linewidth = 0.8) +
  scale_x_continuous(breaks = seq(0, 120, 10)) + 
  labs(x = "Age",
       y = "Probability of death",
       caption = "Mexican females, 2020")



# 05 Save outputs ---------------------------------------------------------

write.csv(x = df_pob_mort_2, file = "data/df_life_table.csv", row.names = FALSE)
ggsave(plot = plt_mort_prob, filename = "figs/plt_mort_prob.png", width = 10, height = 6)

#* Title: Get inflation in medical expenditures in Mexico
#* 
#* Code function: 
#*    This code reads and wrangles the data of the consumer price index specific 
#*    for the  health and personal care area in Mexico from 1990 to 2026
#* 
#* Source of data
#* - https://www.inegi.org.mx/temas/inpc/#tabulados
#* - https://www.inegi.org.mx/app/tabulados/inp/default.aspx?nc=ca57_2018a&idrt=137&opc=t
#* 
#* Creation date: April 02 2026
#* Author: David U. Garibay Treviño

# 01 Initial Setup --------------------------------------------------------

## 01.01 Clean environment ------------------------------------------------
remove(list = ls())

#* Refresh environment memory
gc()

## 01.02 Load libraries ----------------------------------------------------
library(readr)
library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(ggplot2)

# 02 Read data ------------------------------------------------------------

#* Read mexican consumer price index focused in the healthcare sector in Mexico
#* (Indice Nacional de Precios al Consumidor [INPC], clasificación salud y 
#* cuidado personal)
#* 
#* Base index: 4th week of July 2018 = 100
#* Periodicity: Monthly
df_INPC_00 <- read_csv(file = "data-raw/INEGI_INPC_00.csv")




# 03 Wrangle the data -----------------------------------------------------

df_CPI_health <- df_INPC_00 %>% 
  select(date, `health and personal care`) %>% 
  # Convert the date into numeric format
  mutate(date = stringr::str_replace_all(
    string = date,
    pattern = c("Ene " = "01-",
                "Feb " = "02-",
                "Mar " = "03-",
                "Abr " = "04-",
                "May " = "05-",
                "Jun " = "06-",
                "Jul " = "07-",
                "Ago " = "08-",
                "Sep " = "09-",
                "Oct " = "10-",
                "Nov " = "11-",
                "Dic " = "12-"))) %>% 
  tidyr::separate_wider_delim(cols = date, 
                       delim = "-",
                       names = c("month", "year"),
                       cols_remove = FALSE) %>% 
  mutate(month = as.numeric(month),
         year = as.numeric(year)) %>% 
  mutate(date_2 = lubridate::my(date))



ggplot(data = df_CPI_health,
       mapping = aes(x = date_2, 
                     y = `health and personal care`)) + 
  theme_bw() + 
  geom_line() + 
  geom_hline(yintercept = 100,
             linetype = "dotted")
  

# 04 Save the data --------------------------------------------------------
saveRDS(object = df_CPI_health, file = "data/df_CPI_health.rds")

# # Monthly costs
# stage_1 <- 5398/8.38 
# stage_2 <- 4800/8.85
# stage_3 <- 4943/6.38 
# stage_4 <- 5386/5.93 
# 
# mean_cost <- mean(c(stage_1, stage_2, stage_3, stage_4))



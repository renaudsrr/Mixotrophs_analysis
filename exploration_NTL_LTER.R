# Script by SERRE Renaud
# Internship UQAM-GRIL under supervision of BEISNER Beatrix

library(dplyr)
library(openxlsx)
library(tidyverse)
library(gridExtra)
library(corrplot)
library(glue)

library(plotly)


rm(list=ls()) 
graphics.off() 
cat("\014")

setwd("/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/Data/NTL_LTER")

###########################################################################
# Data --------------------------------------------------------------------
###########################################################################

folder_path = "/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/Data/NTL_LTER"
file_list = list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
data_list = lapply(file_list, read.csv)
names(data_list) = gsub("\\.csv$", "", basename(file_list))
print(names(data_list))

info_genus = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/Data/new_csv/info_genus.csv")
all_mixo = info_genus$Genus[!is.na(info_genus$Genus) & info_genus$Final_Nutrition_Strategy == "Mixotroph"]
length(all_mixo)

data_phyto_raw = bind_rows(data_list$phyto_south, data_list$phyto_north)

data_phyto = data_phyto_raw %>%
  select(-c("sampledate","sta","depth_range","division","taxa_name","gald","cells_per_nu","nu_per_ml","cells_per_ml","biovolume_conc","relative_total_biovolume")) %>%
  rename("lake_id" = lakeid,
         "year" = year4,
         "biomass_genus" = biomass_conc) %>%
  mutate(strat = ifelse(genus %in% all_mixo, "Mixotroph", "Autotroph"))

data_phyto2 = data_phyto %>%
  group_by(year,lake_id) %>%
    summarise(
      rich_genus =  n_distinct(genus[biomass_genus>0]), 
      nb_genus_mixo = n_distinct(genus[strat == "Mixotroph" & biomass_genus>0]),
      
      biomass_tot = sum(biomass_genus),
      biomass_mixo = sum(biomass_genus[strat == "Mixotroph"& biomass_genus>0]),
     
      H = {
        x = biomass_genus / sum(biomass_genus)
        x = x[x > 0]
        -sum(x * log(x))},
      `1-D` = {
        x = biomass_genus / sum(biomass_genus)
        sum(x^2)}
       
    ) %>%
    mutate(
    J = ifelse(rich_genus > 1, H / log(rich_genus), NA),
    prev_Mixo = ifelse(biomass_tot > 0, biomass_mixo / biomass_tot * 100, NA))

###########################################################################
# Exploration Data --------------------------------------------------------------------
###########################################################################

choix_div = "J"

fig <- plot_ly()

for (an in unique(data_phyto2$year)) {
  df_year <- data_phyto2 %>%
    filter(year == an, !is.na(!!sym(choix_div)), !is.na(prev_Mixo), !is.na(year)) %>%
    arrange(prev_Mixo)
  
  fig <- fig %>%
    add_trace(x = df_year$prev_Mixo,
              y = df_year$year,
              z = df_year[[choix_div]],
              type = "scatter3d",
              mode = "lines+markers",
              marker = list(size = 3,
                            color="black"),
              line = list(width = 2,
                          color='black'))
}

fig <- fig %>% layout(
  scene = list(
    xaxis = list(title = "prev_mixo (%)"),
    yaxis = list(title = "year"),
    zaxis = list(title = glue("Diversité : {choix_div}"))
  ),
  showlegend = FALSE
)

fig


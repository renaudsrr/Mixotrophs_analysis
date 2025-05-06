# Script by SERRE Renaud
# Internship UQAM-GRIL under supervision of BEISNER Beatrix

rm(list=ls()) 
graphics.off() 
cat("\014")

library(tidyverse)
library(dplyr)
library(corrplot)
library(gridExtra)
library(worrms)
library(purrr)
library(openxlsx)
library(vegan)

#####################################################################################################
# IMPORT DATA ---------------------------------------------------------------------------------------
#####################################################################################################

setwd("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/4_IISD-ELA/csv_traites")
dt1 = read.csv("dt1.csv") # env chimie
dt2 = read.csv("dt2.csv") # infos spatial
dt3 = read.csv("dt3.csv") # biomass plantons
dt4 = read.csv("dt4.csv") # biomass / sp # 
dt7 = read.csv("dt7.csv") # env physics
setwd("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/4_IISD-ELA")
secchi = read.csv("RID_375_field_obs.csv")
temp_and_O2 = read.csv("RID_375_profiles.csv")
temp_surface = read.csv("RID_375_water_surf_temp.csv")

# -----------------------------------------------------------------------------------------------
# Phyto                                                                             
# -----------------------------------------------------------------------------------------------
sp_code = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/4_IISD-ELA/csv_traites/dt5.csv")[,-1]
taxo_info = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/taxo_info.csv")[,-1]

data_phyto = dt4[,c(2,3,7,8:10)] %>% 
  rename(lake_id = monitoring_location_name,
         date_sampling = date_collected) %>%
  mutate(date_sampling = as.Date(date_sampling), format = "%Y-%m-%d") %>%
  left_join(sp_code[,c(1,7)],by="species_code") %>%
  rename(Genus = genus) %>%
  left_join(taxo_info,by="Genus") %>%
  mutate(biovol_um3.mL = (density*volume_cell)/10^6) %>% # calcul biovolume
  dplyr::select(lake_id,date_sampling,species_code,Empire,Kingdom,Phylum,Class,Order,Family,Genus,strat,biovol_um3.mL)
data_phyto = na.omit(data_phyto)

# -----------------------------------------------------------------------------------------------
# Zoo                                                                           
# -----------------------------------------------------------------------------------------------

# A IMPLEMENTER

# -----------------------------------------------------------------------------------------------
# Env                                                                             
# -----------------------------------------------------------------------------------------------

env_chim = dt1 %>%
  dplyr::select(monitoring_location_name,activity_start_date,# activity_depth_height_measure,
         characteristic_name, result_value, result_unit) %>%
  mutate(characteristic_name = replace_na(characteristic_name, "SO")) %>%
  mutate(
     result_unit = str_replace_all(result_unit, "/", "."),
     characteristic_name = paste(characteristic_name, result_unit,sep = "_")) %>% dplyr::select(-c(result_unit)) %>%
  pivot_wider(
    values_from = result_value,
    names_from = characteristic_name,
    values_fn = list(result_value = mean),
    values_fill = list(result_value = NA)) %>%
  rename(lake_id = monitoring_location_name,
         date_sampling = activity_start_date,
         # date_analysis = analysis_start_date,
         # zmax = activity_depth_height_measure,
         pH = PH_,
         SR = SR_,
         COLOUR = COLOUR_,
         DO_mg.L = O2_mg.L
         ) %>%
  mutate(date_sampling = as.Date(date_sampling), format = "%Y-%m-%d") %>%
  group_by(date_sampling,lake_id)%>% # mean by date_sampling, que faire du date_analysis ?????
  summarise(
    mean_ALK_ueq.L = mean(ALK_ueq.L, na.rm = TRUE),
    mean_CA_mg.L = mean(CA_mg.L, na.rm = TRUE),
    mean_CHLA_ug.L = mean(CHLA_ug.L, na.rm = TRUE),
    mean_CL_mg.L = mean(CL_mg.L, na.rm = TRUE),
    mean_COND_uS.cm = mean(COND_uS.cm, na.rm = TRUE),
    mean_DIC_umol.L = mean(DIC_umol.L, na.rm = TRUE),
    mean_DOC_umol.L = mean(DOC_umol.L, na.rm = TRUE),
    mean_FE_mg.L = mean(FE_mg.L, na.rm = TRUE),
    mean_K_mg.L = mean(K_mg.L, na.rm = TRUE),
    mean_MG_mg.L = mean(MG_mg.L, na.rm = TRUE),
    mean_MN_mg.L = mean(MN_mg.L, na.rm = TRUE),
    mean_SO_mg.L = mean(SO_mg.L, na.rm = TRUE),
    mean_NH3_ug.L = mean(NH3_ug.L, na.rm = TRUE),
    mean_NO2_ug.L = mean(NO2_ug.L, na.rm = TRUE),
    mean_NO3_ug.L = mean(NO3_ug.L, na.rm = TRUE),
    mean_PARTC_ug.L = mean(PARTC_ug.L, na.rm = TRUE),
    mean_PARTP_ug.L = mean(PARTP_ug.L, na.rm = TRUE),
    mean_pH = mean(pH, na.rm = TRUE),
    mean_SO4_mg.L = mean(SO4_mg.L, na.rm = TRUE),
    mean_SRSI_mg.L = mean(SRSI_mg.L, na.rm = TRUE),
    mean_TDN_ug.L = mean(TDN_ug.L, na.rm = TRUE),
    mean_TDP_ug.L = mean(TDP_ug.L, na.rm = TRUE),
    mean_PARTN_ug.L = mean(PARTN_ug.L, na.rm = TRUE),
    mean_DO_mg.L = mean(DO_mg.L, na.rm = TRUE),
    mean_PARTFE_ug.L = mean(PARTFE_ug.L, na.rm = TRUE),
    mean_COLOUR = mean(COLOUR, na.rm = TRUE),
    mean_A254_m_m1 = mean(`A254_m^-1`, na.rm = TRUE),
    mean_SR = mean(SR, na.rm = TRUE),
    mean_SRP_ug.L = mean(SRP_ug.L, na.rm = TRUE),
  ) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)))

env_phy = dt7 %>%
  dplyr::select(-c(X,order_lake,mixing_status)) %>%
  rename(lake_id = monitoring_location_name,
         long = longitude,
         lat = latitude)

env = env_chim %>%
  left_join(env_phy,by="lake_id")

#####################################################################################################
# Summary DATA ---------------------------------------------------------------------------------------
#####################################################################################################

# -----------------------------------------------------------------------------------------------
# Phyto                                                                              
# -----------------------------------------------------------------------------------------------

summary_data_phyto = data_phyto %>%
  filter(biovol_um3.mL>0) %>%
  group_by(lake_id,date_sampling) %>%
  summarize(
    biovol_tot_um3.mL= sum(biovol_um3.mL),
    biovol_mixo_um3.mL = sum(biovol_um3.mL[strat == "Mixotroph"]),
    rich_genus = n_distinct(Genus),
    rich_genus_mixo = sum(ifelse(strat=="Mixotroph",1,0)),
    shannon = -sum((biovol_um3.mL/biovol_tot_um3.mL) * log(biovol_um3.mL/biovol_tot_um3.mL) ),
    op_simpson = 1-sum((biovol_um3.mL/biovol_tot_um3.mL)**2)
  ) %>% 
  ungroup() %>%
  mutate(eveness_piel = shannon/log(rich_genus),
         prev_Mixo = biovol_mixo_um3.mL/biovol_tot_um3.mL*100) %>%
  dplyr::select(lake_id,date_sampling,everything())

# -----------------------------------------------------------------------------------------------
# Phyto corrigé                                                                        
# -----------------------------------------------------------------------------------------------

summary_data_phyto_corr = data_phyto %>%
  filter(biovol_um3.mL>0) %>%
  filter(strat=="Autotroph") %>%
  
  group_by(lake_id,date_sampling) %>%
  summarize(
    biovol_tot_um3.mL= sum(biovol_um3.mL),
    
    rich_genus_corr = n_distinct(Genus),
    shannon_corr = -sum((biovol_um3.mL/biovol_tot_um3.mL) * log(biovol_um3.mL/biovol_tot_um3.mL) ),
    op_simpson_corr = 1-sum((biovol_um3.mL/biovol_tot_um3.mL)**2)
  ) %>% 
  ungroup() %>%
  mutate(eveness_piel_corr = shannon_corr/log(rich_genus_corr)) %>%
  dplyr::select(lake_id,date_sampling,everything())

# -----------------------------------------------------------------------------------------------
# Zoo                                                                             
# -----------------------------------------------------------------------------------------------

# A IMPLEMENTER

#####################################################################################################
# DIV BETA  ---------------------------------------------------------------------------------------
#####################################################################################################

# A IMPLEMENTER

#####################################################################################################
# Export CSV ---------------------------------------------------------------------------------------
#####################################################################################################

data_python = env %>%
  left_join(summary_data_phyto, by = c("lake_id", "date_sampling")) %>%
  left_join(summary_data_phyto_corr[,-3],by=c("lake_id", "date_sampling")) %>%
  filter(!is.na(shannon)) %>% # correspondance entre env et data_phyto
  # left_join(summary_data_zoo[], by = c("lake_id", "date_sampling")) %>%
  mutate(doy = yday(date_sampling),
         yr = year(date_sampling)) %>%
  dplyr::select(-c("date_sampling","biovol_tot_um3.mL","biovol_mixo_um3.mL","rich_genus_mixo"))
data_python = data_python[,-1] # car "> Adding missing grouping variables: `date_sampling`

data_python = data_python[!is.na(data_python$shannon), ]
data_python = data_python[!is.na(data_python$eveness_piel), ]

write.csv(data_python,"/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python/df_IIDS_ELA.csv",row.names=F)

#####################################################################################################
# Analyses ---------------------------------------------------------------------------------------
#####################################################################################################
# 
# # --- All VA en fonction de prev_mixo ---------------------------------------------------------------------------
# numeric_vars = data_python %>%
#   dplyr::select(where(is.numeric)) %>%
#   dplyr::select(-prev_Mixo) %>%
#   colnames()
# 
# for (var in numeric_vars) {
#   gg = ggplot(data_python, aes_string(x = "prev_Mixo", y = var)) +
#     geom_point(alpha = 0.5) +
#     geom_smooth(method = "gam", se = T, color = "red") + # ajuste un gam sur chaque var
#     labs(x = "prev_Mixo",
#          y = var) +
#     theme_minimal()
# 
#   print(gg)
# }
# 
# # --- All VA en fonction de doy ---------------------------------------------------------------------------
# numeric_vars = data_python %>%
#   dplyr::select(where(is.numeric)) %>%
#   dplyr::select(-doy) %>%
#   colnames()
# 
# for (var in numeric_vars) {
#   gg = ggplot(data_python, aes_string(x = "doy", y = var)) +
#     geom_point(alpha = 0.5) +
#     geom_smooth(method = "gam", se = T, color = "red") + # ajuste un gam sur chaque var
#     labs(x = "doy",
#          y = var) +
#     theme_minimal()
#   
#   print(gg)
# }
# 
# # --- All VA en fonction de doy_cumul ---------------------------------------------------------------------------
# numeric_vars = data_python %>%
#   dplyr::select(where(is.numeric)) %>%
#   dplyr::select(-doy_cumul) %>%
#   colnames()
# 
# for (var in numeric_vars) {
#   gg = ggplot(data_python, aes_string(x = "doy_cumul", y = var)) +
#     geom_point(alpha = 0.5) +
#     geom_smooth(method = "gam", se = T, color = "red") + # ajuste un gam sur chaque var
#     labs(x = "doy_cumul",
#          y = var) +
#     theme_minimal()
#   
#   print(gg)
# }



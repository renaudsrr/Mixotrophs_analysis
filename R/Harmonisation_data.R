# Script by SERRE Renaud
# Internship UQAM-GRIL under supervision of BEISNER Beatrix

rm(list=ls()) 
graphics.off() 
cat("\014")

library(dplyr)
library(openxlsx)
library(tidyverse)
library(gridExtra)
library(grid) 
library(lubridate)

setwd("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python")

data_phyton_CARBBAS = read.csv("df_CARBBAS.csv")
data_phyton_LP_NLA = read.csv("df_LP_NLA.csv")
data_python_IISD_ELA = read.csv("df_IIDS_ELA.csv")

colnames(data_phyton_CARBBAS)
colnames(data_phyton_LP_NLA)
colnames(data_python_IISD_ELA)

# --- CARBBAS ---------------------------------------------------------------------------
data_phyton_CARBBAS_harmz = data_phyton_CARBBAS %>%
  rename(DOC_mg.L = doc.mgl,
         Chla_ug.L = chl.ugl,
         sech_m = sech.m,
         temp_mean = h2otemp,
         pH = ph,
         DO_mg.L = do.mgl,
         TP_ug.L = tp.ugl,
         area_m2 = ctch.m2,
         depth = zmax) %>%
  mutate(TN_ug.L = tn.mgl*1000) %>%
  dplyr::select(lake_id,yr,doy,lat,long,area_m2,depth,color,sech_m,temp_mean,pH,Chla_ug.L,DO_mg.L,DOC_mg.L,TP_ug.L,TN_ug.L,
         prev_Mixo,biovol_tot_zoo_um3.mL,biovol_cladocera_um3.mL,biovol_copepoda_um3.mL,biovol_other_um3.mL,
         rich_genus,shannon,op_simpson,eveness_piel)

# --- LP-NLA rapide ---------------------------------------------------------------------------
data_phyton_LP_NLA_harmz = data_phyton_LP_NLA %>%
  rename(DOC_mg.L = DOC, 
         color = Colour,
         Chla_ug.L = Chla, 
         sech_m = secchi_depth,
         lat = latitude,
         long = longitude,
         area_m2 = area,
         depth = lake_depth,
         SO_mg.L = Sodium,
         MG_mg.L = Magnesium,
         CL_mg.L = Chloride,
         K_mg.L = Potassium,
         TP_ug.L = TP
         ) %>%
  mutate(temp_mean = rowMeans(cbind(Temperature_up, Temperature_bottom),na.rm=T),
         pH = rowMeans(cbind(Temperature_up, Temperature_bottom),na.rm = T),
         DO_mg.L = rowMeans(cbind(Temperature_up, Temperature_bottom),na.rm = T),
         TN_ug.L = TN*1000, # car en mg/L
         COND_uS.cm = Conductivity*1000 # car en mS/cm
         ) %>%
  dplyr::select(lake_id,lat,long,area_m2,depth,color,sech_m,temp_mean,pH,Chla_ug.L,DO_mg.L,DOC_mg.L,TP_ug.L,TN_ug.L,SO_mg.L,MG_mg.L,CL_mg.L,K_mg.L,COND_uS.cm, 
         biovol_tot_zoo_um3.mL,biovol_cladocera_um3.mL,biovol_copepoda_um3.mL,biovol_other_um3.mL,
         rich_genus,shannon,op_simpson,eveness_piel,log_LCBD,
         rich_genus_corr,shannon_corr,op_simpson_corr,eveness_piel_corr,log_LCBD_corr,
         prev_Mixo
         )

# --- IISD-ELA ---------------------------------------------------------------------------
data_phyton_IISD_ELA_harmz = data_python_IISD_ELA %>%
  rename(
    color = mean_COLOUR,
    Chla_ug.L = mean_CHLA_ug.L,
    pH = mean_pH,
    DO_mg.L = mean_DO_mg.L,
    TP_ug.L = mean_TDP_ug.L,
    TN_ug.L = mean_TDN_ug.L,
    area_m2 = area_surface,
    depth = depth_max,
    SO_mg.L = mean_SO_mg.L,
    MG_mg.L = mean_MG_mg.L,
    CL_mg.L = mean_CL_mg.L,
    K_mg.L = mean_K_mg.L,
    COND_uS.cm = mean_COND_uS.cm
  ) %>%
  mutate(
    DOC_mg.L = mean_DOC_umol.L * 0.012, # car en umol/L
    # PROVISOIRE-----------------
    biovol_tot_zoo_um3.mL = 0,
    biovol_cladocera_um3.mL = 0,
    biovol_copepoda_um3.mL = 0,
    biovol_other_um3.mL = 0,
    sech_m = 0,
    temp_mean = 0
    # PROVISOIRE-----------------
  ) %>%
  dplyr::select(lake_id,yr,doy,lat,long,area_m2,depth,color,sech_m,temp_mean,pH,Chla_ug.L,DO_mg.L,DOC_mg.L,TP_ug.L,TN_ug.L,SO_mg.L,MG_mg.L,CL_mg.L,K_mg.L,COND_uS.cm,
         biovol_tot_zoo_um3.mL,biovol_cladocera_um3.mL,biovol_copepoda_um3.mL,biovol_other_um3.mL,
         rich_genus,shannon,op_simpson,eveness_piel,
         rich_genus_corr,shannon_corr,op_simpson_corr,eveness_piel_corr,
         prev_Mixo)


# write csv ---------------------------------------------------------------
setwd("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python")

write.csv(data_phyton_CARBBAS_harmz,"1_harmz_CARBBAS.csv",row.names=F)
write.csv(data_phyton_LP_NLA_harmz,"2_harmz_LP_NLA.csv",row.names=F)
write.csv(data_phyton_IISD_ELA_harmz,"3_harmz_IISD_ELA.csv",row.names=F)





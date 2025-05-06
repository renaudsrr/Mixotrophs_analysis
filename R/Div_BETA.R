# Script by SERRE Renaud
# Internship UQAM-GRIL under supervision of BEISNER Beatrix

rm(list=ls()) 
graphics.off() 
cat("\014")

library(dplyr)
library(ggplot2)
library(tidyr)
library(vegan)
library(lubridate)

#####################################################################################################
# DIVERSITE BETA LP-NLA -----------------------------------------------------------------------------
#####################################################################################################

data_phyto = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/2_LP_NLA/HarmonizedPhyto_LakePulseNLA2017_19092022.csv")

# --- Ytr ---------------------------------------------------------------------------
Y = data_phyto %>%
  rename(Genus = Target_taxon) %>%
  dplyr::select(lake_id,Genus,biovol) %>%
  pivot_wider(names_from = Genus,values_from = biovol,values_fill = 0)
lake_ids = Y$lake_id
Y = as.matrix(Y[,-1])  # n sites x p genus
rownames(Y) = lake_ids
Ytr = decostand(Y, method = "hellinger") 

# --- S ---------------------------------------------------------------------------
S = matrix(nrow = nrow(Ytr),ncol=ncol(Ytr),dimnames=list(rownames(Ytr),colnames(Ytr)))

for (col in 1:ncol(S)){
  for (row in 1:nrow(S)){
    mean_yj = mean(Ytr[,col])
    yij = Ytr[row,col]
    sij = (yij-mean_yj)^2
    S[row,col] = sij
  }
}
# Sbis = (scale(Ytr, center = T, scale = F))^2

# --- Beta div ---------------------------------------------------------------------------
SS_Total = sum(S)
BD_Total = SS_Total / (nrow(S)-1)

# Contrib de chaque genus à la BD_Total
SSj = colSums(S)
df_SCBD = as.data.frame(SSj / SS_Total)
colnames(df_SCBD) = "SCBD"
df_SCBD$Genus = as.factor(rownames(df_SCBD))
rownames(df_SCBD) = NULL

# COntrib de chaque lakes à la BD_Total
SSi = rowSums(S)
df_LCBD = as.data.frame(SSi / SS_Total)
colnames(df_LCBD) = "LCBD"
df_LCBD$lake_id = as.factor(rownames(df_LCBD))
rownames(df_LCBD) = NULL

# -----------------------------------------------------------------------------------------------
#  Plot et analyse                                                                             
# -----------------------------------------------------------------------------------------------

# --- importance genus ---------------------------------------------------------------------------
thresh = 0.01
df_thresh = df_SCBD[which(df_SCBD$SCBD>thresh),]
barplot(names.arg=droplevels(df_thresh$Genus),height=df_thresh$SCBD,las=2)

# --- importance lakes ---------------------------------------------------------------------------  
thresh = 8e-04
df_thresh = df_LCBD[which(df_LCBD$LCBD>thresh),]
barplot(names.arg=droplevels(df_thresh$lake_id),height=df_thresh$LCBD,las=2)

# --- Test de permutation ---------------------------------------------------------------------------
# # H0 : les espèces sont réparties aléatoirement entre tout les sites, il n'y a pas d'organisation écologique
# 
# nperm = 999 
# 
# LCBD_perm = matrix(NA, nrow = nrow(Ytr), ncol = nperm)
# 
# for (i in 1:nperm) {
#   Y_perm = apply(Ytr, 2, sample)
#   
#   S_perm = (scale(Y_perm, center = TRUE, scale = FALSE))^2
#   
#   SS_Total_perm = sum(S_perm)
#   SSi_perm = rowSums(S_perm)
#   
#   LCBD_perm[, i] = SSi_perm / SS_Total_perm
# }
# 
# # Calcul des p-values empiriques
# LCBD_obs = df_LCBD$LCBD  # LCBD observés
# 
# pval_LCBD = sapply(1:length(LCBD_obs), function(i) {
#   (sum(LCBD_perm[i, ] >= LCBD_obs[i]) + 1) / (nperm + 1)
# })
# 
# df_LCBD$pval = pval_LCBD
# which(df_LCBD$pval<0.05)
# Aucun lac ne se détache -> pas de structure écologique forte détectée

# --- Regression ---------------------------------------------------------------------------
data_python = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python/2_harmz_LP_NLA.csv") %>%
  dplyr::select(-c("LCBD"))

df_beta = data_python %>%
  left_join(df_LCBD[,c("lake_id","LCBD")],by="lake_id")

plot(df_beta$prev_Mixo,df_beta$LCBD)

# --- Div LCBD ~ Predictors ---------------------------------------------------------------------------
# numeric_vars = df_beta %>%
#   dplyr::select(where(is.numeric))
# 
# for (var in names(numeric_vars)) {
#   gg = ggplot(df_beta, aes(x = .data[[var]], y = LCBD)) +
#     geom_point(alpha = 0.5) +
#     geom_smooth(method = "gam", se = T, color = "red") + # ajuste un gam sur chaque var
#     labs(x = var,
#          y = "LCBD") +
#     theme_minimal()
# 
#   print(gg)
# }

# --- Diag Venn ---------------------------------------------------------------------------
data_python = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python/2_harmz_LP_NLA.csv",row.names = 1)
VA_space = data_python[,c("long","lat")]
VA_env = data_python %>% 
  dplyr::select(-c("long","lat")) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))
lacs_communs = intersect(rownames(VA_env), rownames(Ytr))
SP = Ytr[lacs_communs, ]
res_varpart = varpart(SP, VA_space, VA_env)
plot(res_varpart, digits = 2)



#####################################################################################################
# DIVERSITE BETA ELA --------------------------------------------------------------------------------
#####################################################################################################

dt4 = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/4_IISD-ELA/csv_traites/dt4.csv")
sp_code = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/4_IISD-ELA/csv_traites/dt5.csv")[,-1]
taxo_info = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/taxo_info.csv")[,-1]

data_phyto2 = dt4[,c(2,3,7,8:10)] %>% 
  rename(lake_id = monitoring_location_name,
         date_sampling = date_collected) %>%
  mutate(date_sampling = as.Date(date_sampling), format = "%Y-%m-%d") %>%
  left_join(sp_code[,c(1,7)],by="species_code") %>%
  rename(Genus = genus) %>%
  left_join(taxo_info,by="Genus") %>%
  mutate(biovol_um3.mL = (density*volume_cell)/10^6) %>% # calcul biovolume
  dplyr::select(lake_id,date_sampling,species_code,Empire,Kingdom,Phylum,Class,Order,Family,Genus,strat,biovol_um3.mL)
data_phyto2 = na.omit(data_phyto2) 

data_phyto2 = data_phyto2 %>%
  mutate(doy = yday(date_sampling)) %>%
  mutate(lake_date = paste(lake_id, doy, sep = "_")) %>%
  dplyr::select(lake_date,Genus,biovol_um3.mL)

data_phyto2 = data_phyto2 %>%
  group_by(lake_date, Genus) %>%
  summarise(biovol_um3.mL = sum(biovol_um3.mL, na.rm = TRUE)) %>%
  ungroup()

Y = data_phyto2 %>%
  pivot_wider(names_from = Genus,values_from = biovol_um3.mL,values_fill = 0)

lake_ids = Y$lake_date
Y = as.matrix(Y[,-1])  # n sites x p genus
rownames(Y) = lake_ids
Ytr = decostand(Y, method = "hellinger") 

# --- S ---------------------------------------------------------------------------
S = matrix(nrow = nrow(Ytr),ncol=ncol(Ytr),dimnames=list(rownames(Ytr),colnames(Ytr)))

# for (col in 1:ncol(S)){
#   for (row in 1:nrow(S)){
#     mean_yj = mean(Ytr[,col])
#     yij = Ytr[row,col]
#     sij = (yij-mean_yj)^2
#     S[row,col] = sij
#   }
# }
S = (scale(Ytr, center = T, scale = F))^2

# --- Beta div ---------------------------------------------------------------------------
SS_Total = sum(S)
BD_Total = SS_Total / (nrow(S)-1)

# Contrib de chaque genus à la BD_Total
SSj = colSums(S)
df_SCBD = as.data.frame(SSj / SS_Total)
colnames(df_SCBD) = "SCBD"
df_SCBD$Genus = as.factor(rownames(df_SCBD))
rownames(df_SCBD) = NULL

# COntrib de chaque lakes à la BD_Total
SSi = rowSums(S)
df_LCBD = as.data.frame(SSi / SS_Total)
colnames(df_LCBD) = "LCBD"
df_LCBD$lake_date = as.factor(rownames(df_LCBD))
rownames(df_LCBD) = NULL

# -----------------------------------------------------------------------------------------------
#  Plot et analyse                                                                             
# -----------------------------------------------------------------------------------------------

# --- importance genus ---------------------------------------------------------------------------
thresh = 0.01
df_thresh = df_SCBD[which(df_SCBD$SCBD>thresh),]
barplot(names.arg=droplevels(df_thresh$Genus),height=df_thresh$SCBD,las=2)

# --- importance lakes ---------------------------------------------------------------------------  
thresh = 8e-04
df_thresh = df_LCBD[which(df_LCBD$LCBD>thresh),]
barplot(names.arg=droplevels(df_thresh$lake_date),height=df_thresh$LCBD,las=2)

# --- Test de permutation ---------------------------------------------------------------------------
# H0 : les espèces sont réparties aléatoirement entre tout les sites, il n'y a pas d'organisation écologique
# 
# nperm = 999 
# 
# LCBD_perm = matrix(NA, nrow = nrow(Ytr), ncol = nperm)
# 
# for (i in 1:nperm) {
#   Y_perm = apply(Ytr, 2, sample)
#   
#   S_perm = (scale(Y_perm, center = TRUE, scale = FALSE))^2
#   
#   SS_Total_perm = sum(S_perm)
#   SSi_perm = rowSums(S_perm)
#   
#   LCBD_perm[, i] = SSi_perm / SS_Total_perm
# }

# Calcul des p-values empiriques
# LCBD_obs = df_LCBD$LCBD  # LCBD observés
# 
# pval_LCBD = sapply(1:length(LCBD_obs), function(i) {
#   (sum(LCBD_perm[i, ] >= LCBD_obs[i]) + 1) / (nperm + 1)
# })
# 
# df_LCBD$pval = pval_LCBD
# which(df_LCBD$pval<0.05)



# --- Regression ---------------------------------------------------------------------------
data_python = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python/3_harmz_IISD_ELA.csv")

df_beta = data_python %>%
  mutate(lake_date = paste(lake_id, doy, sep = "_")) %>%
  left_join(df_LCBD[,c("lake_date","LCBD")],by="lake_date")

plot(df_beta$prev_Mixo,df_beta$LCBD)

# --- Div LCBD ~ Predictors ---------------------------------------------------------------------------
numeric_vars = df_beta %>%
  dplyr::select(where(is.numeric))

for (var in names(numeric_vars)) {
  gg = ggplot(df_beta, aes(x = .data[[var]], y = LCBD)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "gam", se = T, color = "red") + # ajuste un gam sur chaque var
    labs(x = var,
         y = "LCBD") +
    theme_minimal()

  print(gg)
}

# --- Diag Venn ---------------------------------------------------------------------------
data_python = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python/2_harmz_LP_NLA.csv",row.names = 1)
VA_space = data_python[,c("long","lat")]
VA_env = data_python %>% 
  dplyr::select(-c("long","lat")) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))
lacs_communs = intersect(rownames(VA_env), rownames(Ytr))
SP = Ytr[lacs_communs, ]
# res_varpart = varpart(SP, VA_space, VA_env)
# plot(res_varpart, digits = 2)






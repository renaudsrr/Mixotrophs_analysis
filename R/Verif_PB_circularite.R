# Script by SERRE Renaud
# Internship UQAM-GRIL under supervision of BEISNER Beatrix

rm(list=ls()) 
graphics.off() 
cat("\014")

library(dplyr)
library(ggplot2)
library(tidyr)
library(mgcv) # GAM
library(MASS)
library(DHARMa)
library(rcompanion)
library(vegan)

#####################################################################################################
# DATA ---------------------------------------------------------------------------------------
#####################################################################################################

data_phyto = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/2_LP_NLA/HarmonizedPhyto_LakePulseNLA2017_19092022.csv")
taxo_info = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/taxo_info.csv")[,-1]
all_mixo = taxo_info$Genus[!is.na(taxo_info$Genus) & taxo_info$strat == "Mixotroph"]

summary_data_phyto = data_phyto %>%
  rename("genus" = "Target_taxon",
         "biovol_tot" = "total.biov") %>%
  dplyr::select(-c("classic.group"))  %>%
  mutate(strat = ifelse(genus %in% all_mixo, "Mixotroph", "Autotroph")) %>%
  filter(biovol > 0) %>%
  group_by(lake_id)%>%
  summarize(
    Survey = first(Survey),
    biovol_tot_um3.mL = first(biovol_tot),
    biovol_mixo_um3.mL = sum(biovol[strat == "Mixotroph"]),
    rich_genus = n_distinct(genus),
    rich_genus_mixo = sum(ifelse(strat=="Mixotroph",1,0)),
    shannon = -sum((biovol/biovol_tot_um3.mL) * log(biovol/biovol_tot_um3.mL) ),
    op_simpson = 1-sum((biovol/biovol_tot_um3.mL)**2),
  ) %>%
  ungroup() %>%
  mutate(eveness_piel = shannon/log(rich_genus),
         prev_Mixo = biovol_mixo_um3.mL/biovol_tot_um3.mL*100) 

#####################################################################################################
# VERIFICATION RELATION = ARTEFACT STRUCTUREL ? -----------------------------------------------------
#####################################################################################################

# -----------------------------------------------------------------------------------------------
# Permutation                                                                             
# -----------------------------------------------------------------------------------------------

set.seed(123)
nperm = 1000

summary_data_phyto_obs = summary_data_phyto %>%
  dplyr::select(shannon, prev_Mixo)
mod_obs = gam(shannon ~ s(prev_Mixo), data = summary_data_phyto_obs)
r2_obs = summary(mod_obs)$r.sq

i = 1
all_r2_perm = numeric(nperm)

perm = F
if (perm == T){
  while (i<=nperm) {
    print(paste("perm :",i))
    
    summary_data_phyto_perm = summary_data_phyto_obs %>%
      mutate(prev_Mixo = sample(prev_Mixo))
    
    mod = gam(shannon ~ s(prev_Mixo), data = summary_data_phyto_perm)
    r2 = summary(mod)$r.sq
    all_r2_perm[i] = r2
    
    i = i + 1 
    
  }
  abline(v=r2_obs,col="red",lwd=2)
  hist(all_r2_perm,breaks=200)
  
  pval = sum(all_r2_perm>=r2_obs)/nperm
}

# -----------------------------------------------------------------------------------------------
# Calcul div sans genres mixotrophes                                                                             
# -----------------------------------------------------------------------------------------------

summary_data_phyto_corr = data_phyto %>%
  rename("genus" = "Target_taxon",
         "biovol_tot" = "total.biov") %>%
  dplyr::select(-c("classic.group","biovol_tot"))  %>%
  
  mutate(strat = ifelse(genus %in% all_mixo, "Mixotroph", "Autotroph")) %>%
  filter(biovol > 0) %>%
  filter(strat=="Autotroph") %>%
  
  group_by(lake_id) %>%
  summarize(
    biovol_tot_um3.mL = sum(biovol),
    
    rich_genus_corr = n_distinct(genus),
    shannon_corr = -sum((biovol/biovol_tot_um3.mL) * log(biovol/biovol_tot_um3.mL) ),
    op_simpson_corr = 1-sum((biovol/biovol_tot_um3.mL)**2),
    
  ) %>%
  ungroup() %>%
  mutate(eveness_piel_corr = shannon_corr/log(rich_genus_corr)) %>%
  left_join(summary_data_phyto[,c(1,10)],by=c("lake_id"))

plot(x=summary_data_phyto_corr$prev_Mixo, y = summary_data_phyto_corr$rich_genus_corr)
plot(x=summary_data_phyto_corr$prev_Mixo, y = summary_data_phyto_corr$shannon_corr)
plot(x=summary_data_phyto_corr$prev_Mixo, y = summary_data_phyto_corr$op_simpson_corr)
plot(x=summary_data_phyto_corr$prev_Mixo, y = summary_data_phyto_corr$eveness_piel_corr)


ggplot(summary_data_phyto_corr,aes(x=prev_Mixo,y=rich_genus_corr)) +
  geom_point()+
  geom_smooth(method = "lm",col="red") 
ggplot(summary_data_phyto_corr,aes(x=prev_Mixo,y=shannon_corr)) +
  geom_point()+
  geom_smooth(method = "lm",col="red") 
ggplot(summary_data_phyto_corr,aes(x=prev_Mixo,y=op_simpson_corr)) +
  geom_point()+
  geom_smooth(method = "lm",col="red") 
ggplot(summary_data_phyto_corr,aes(x=prev_Mixo,y=eveness_piel_corr)) +
  geom_point()+
  geom_smooth(method = "lm",col="red") 

#####################################################################################################
# DIV BETA corr ---------------------------------------------------------------------------------------
#####################################################################################################

# --- Ytr ---------------------------------------------------------------------------
Y = data_phyto %>%
  rename(Genus = Target_taxon) %>%
  mutate(strat = ifelse(Genus %in% all_mixo, "Mixotroph", "Autotroph")) %>%
  filter(strat=="Autotroph") %>%
  dplyr::select(lake_id,Genus,biovol) %>%
  pivot_wider(names_from = Genus,values_from = biovol,values_fill = 0)
lake_ids = Y$lake_id
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
df_LCBD$lake_id = as.factor(rownames(df_LCBD))
rownames(df_LCBD) = NULL


data_python = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python/2_harmz_LP_NLA.csv")

df_beta = data_python %>%
  left_join(df_LCBD[,c("lake_id","LCBD")],by="lake_id")

ggplot(df_beta,aes(x=prev_Mixo,y=LCBD)) +
  geom_point()+
  geom_smooth(method = "gam",col="red") 




# IMPORTANCE DE LA RICHESSE MAINTENANTT !!
#####################################################################################################
# Modele sur la richesse spé ---------------------------------------------------------------------------------------
#####################################################################################################

# --- BinomNeg ---------------------------------------------------------------------------

data_python = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python/df_LP_NLA.csv")
colSums(is.na(data_python))
data_python = data_python %>%
  dplyr::select(-c(Temperature_bottom,Thermocline_m,pH_bottom,Epilimnion_bot,Hypolimnion_top,biovol_tot_zoo_um3.mL,biovol_cladocera_um3.mL,biovol_copepoda_um3.mL,biovol_other_um3.mL))
colSums(is.na(data_python))

data_python = data_python %>%
  left_join(df_LCBD[,c("lake_id","LCBD")],by="lake_id") %>%
  left_join(summary_data_phyto_corr[,c(1,3:6)],by="lake_id")

nrow(data_python)
data_python = na.omit(data_python)
nrow(data_python)









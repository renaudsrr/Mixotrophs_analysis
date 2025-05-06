# Script by SERRE Renaud
# Internship UQAM-GRIL under supervision of BEISNER Beatrix

rm(list=ls()) 
graphics.off() 
cat("\014")

library(dplyr)
library(ggplot2)
library(tidyr)
library(corrplot)
library(gridExtra)
library(worrms)
library(purrr)
library(openxlsx)
library(vegan)
library(tibble)

#####################################################################################################
# IMPORT DATA ---------------------------------------------------------------------------------------
#####################################################################################################

# -----------------------------------------------------------------------------------------------
# Phyto                                                                             
# -----------------------------------------------------------------------------------------------

data_phyto = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/2_LP_NLA/HarmonizedPhyto_LakePulseNLA2017_19092022.csv")

# -----------------------------------------------------------------------------------------------
# Zoo                                                                           
# -----------------------------------------------------------------------------------------------

data_zoo = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/2_LP_NLA/ZooLakePulseRaw_ALLfinal_grouping2017_2018_2019.csv",sep=";") %>%
  rename("lake_id"="Lake_ID")
colnames(data_zoo) = sub("\\..*", "", colnames(data_zoo))
colnames(data_zoo) = make.names(colnames(data_zoo), unique = TRUE)
data_zoo = data_zoo %>%
  mutate(
    Daphnia = Daphnia + Daphnia.1 + Daphnia.2 + Daphnia.3 + Daphnia.4 + Daphnia.5 + Daphnia.6 + Daphnia.7 + Daphnia.8 + Daphnia.9 + Daphnia.10 + Daphnia.11 + Daphnia.12 + Daphnia.13,
    Acanthocyclops = Acanthocyclops + Acanthocyclops.1,
    Eucyclops = Eucyclops + Eucyclops.1 + Eucyclops.2,
    Mesocyclops = Mesocyclops + Mesocyclops.1,
    Epischura = Epischura + Epischura.1 + Epischura.2 + Epischura.3,
    Aglaodiaptomus = Aglaodiaptomus + Aglaodiaptomus.1 + Aglaodiaptomus.2 + Aglaodiaptomus.3,
    Leptodiaptomus = Leptodiaptomus + Leptodiaptomus.1 + Leptodiaptomus.2 + Leptodiaptomus.3 + Leptodiaptomus.4 + Leptodiaptomus.5 + Leptodiaptomus.6,
    Skistodiaptomus = Skistodiaptomus + Skistodiaptomus.1 + Skistodiaptomus.2 + Skistodiaptomus.3 + Skistodiaptomus.4,
    Alona = Alona + Alona.1 + Alona.2 + Alona.3 + Alona.4 + Alona.5,
    Alonella = Alonella + Alonella.1,
    Ceriodaphnia = Ceriodaphnia + Ceriodaphnia.1 + Ceriodaphnia.2,
    Diaphanosoma = Diaphanosoma + Diaphanosoma.1 + Diaphanosoma.2,
    Eubosmina = Eubosmina + Eubosmina.1 + Eubosmina.2 + Eubosmina.3,
    Hesperodiaptomus = Hesperodiaptomus + Hesperodiaptomus.1 + Hesperodiaptomus.2 + Hesperodiaptomus.3 + Hesperodiaptomus.4 + Hesperodiaptomus.5,
    Onychodiaptomus = Onychodiaptomus + Onychodiaptomus.1,
    Cyclops = Cyclops + Cyclops.1,
    Diacyclops = Diacyclops + Diacyclops.1
  ) %>%
  dplyr::select(-c(
    Diacyclops.1,
    Cyclops.1,
    Daphnia.1, Daphnia.2, Daphnia.3, Daphnia.4, Daphnia.5, Daphnia.6, Daphnia.7, Daphnia.8, Daphnia.9, Daphnia.10, Daphnia.11, Daphnia.12, Daphnia.13,
    Acanthocyclops.1,
    Eucyclops.1, Eucyclops.2,
    Mesocyclops.1,
    Epischura.1, Epischura.2, Epischura.3,
    Aglaodiaptomus.1, Aglaodiaptomus.2, Aglaodiaptomus.3,
    Leptodiaptomus.1, Leptodiaptomus.2, Leptodiaptomus.3, Leptodiaptomus.4, Leptodiaptomus.5, Leptodiaptomus.6,
    Skistodiaptomus.1, Skistodiaptomus.2, Skistodiaptomus.3, Skistodiaptomus.4,
    Alona.1, Alona.2, Alona.3, Alona.4, Alona.5,
    Alonella.1,
    Ceriodaphnia.1, Ceriodaphnia.2,
    Diaphanosoma.1, Diaphanosoma.2,
    Eubosmina.1, Eubosmina.2, Eubosmina.3,
    Hesperodiaptomus.1, Hesperodiaptomus.2, Hesperodiaptomus.3, Hesperodiaptomus.4, Hesperodiaptomus.5,
    Onychodiaptomus.1
  ))

# -----------------------------------------------------------------------------------------------
# Env                                                                             
# -----------------------------------------------------------------------------------------------
env = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/2_LP_NLA/HarmonizedPredictors_LakePulseNLA2017_21062022.csv") %>%
  rename("lake_id"="Lake_ID") %>% dplyr::select(Survey,everything())

# -----------------------------------------------------------------------------------------------
# Infos genus                                                                             
# -----------------------------------------------------------------------------------------------
taxo_info = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/taxo_info.csv")[,-1]
all_mixo = taxo_info$Genus[!is.na(taxo_info$Genus) & taxo_info$strat == "Mixotroph"]

#####################################################################################################
# Summary DATA ---------------------------------------------------------------------------------------
#####################################################################################################

# -----------------------------------------------------------------------------------------------
# Phyto                                                                             
# -----------------------------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------------------------
# Phyto corrigé (Autotrophs)                                                                             
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
  mutate(eveness_piel_corr = shannon_corr/log(rich_genus_corr))

# -----------------------------------------------------------------------------------------------
# Zoo                                                                             
# -----------------------------------------------------------------------------------------------
# Classification en 3 groupes : Cladocera, Copepoda et OTHER
cladocera = c("Daphnia", "Acroperus", "Alona", "Alonella", "Bosmina", "Camptocercus", "Ceriodaphnia", "Chydorus",
              "Diaphanosoma", "Eubosmina", "Holopedium", "Macrothrix", "Picripleuroxus", "Pleuroxus", "Polyphemus",
              "Sida", "Simocephalus", "Dunhevedia", "Graptoleberis", "Ilyocryptus", "Latona", "Moina", "Eurycercus",
              "Scapholeberis")

copepoda = c("Ergasilus", "harpacticoid", "Acanthocyclops", "Cyclops", "Diacyclops", "Eucyclops", "Mesocyclops",
             "Microcyclops", "Orthocyclops", "Tropocyclops", "cyclopoid", "Aglaodiaptomus", "Epischura",
             "Leptodiaptomus", "Limnocalanus", "Onychodiaptomus", "Senecella", "Skistodiaptomus", "calanoid",
             "Diacyclops", "Macrocyclops", "Acanthodiaptomus", "Hesperodiaptomus", "Cyclops", "Heterocope")

other = c("ostracod")

H = diversity(as.matrix(data_zoo[,-1]), index="shannon")
one_minus_D = diversity(as.matrix(data_zoo[,-1]), index="simpson") # return 1-D

summary_data_zoo = data_zoo %>%
  rowwise() %>%
  mutate(
    biovol_cladocera_um3.mL = sum(c_across(all_of(cladocera)), na.rm = TRUE),
    biovol_copepoda_um3.mL = sum(c_across(all_of(copepoda)), na.rm = TRUE),
    biovol_other_um3.mL = sum(c_across(all_of(other)), na.rm = TRUE),
    biovol_tot_zoo_um3.mL = sum(c_across(-lake_id), na.rm = TRUE),
    rich_genus = sum(c_across(!any_of("lake_id")) > 0)
  ) %>%
  ungroup() %>%
  mutate(shannon = H,
         op_simpson = one_minus_D, 
         eveness_piel = H/log(rich_genus)) %>%
  left_join(summary_data_phyto[,c(1,10)],by="lake_id") %>%
  dplyr::select(lake_id, biovol_tot_zoo_um3.mL, 
         biovol_cladocera_um3.mL, biovol_copepoda_um3.mL, biovol_other_um3.mL, 
         rich_genus, shannon, op_simpson, eveness_piel, prev_Mixo)


#####################################################################################################
# DIV BETA ---------------------------------------------------------------------------------------
#####################################################################################################

# -----------------------------------------------------------------------------------------------
# Phyto                                                                             
# -----------------------------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------------------------
# Phyto corr                                                                             
# -----------------------------------------------------------------------------------------------

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
df_SCBD_corr = as.data.frame(SSj / SS_Total)
colnames(df_SCBD_corr) = "SCBD_corr"
df_SCBD_corr$Genus = as.factor(rownames(df_SCBD_corr))
rownames(df_SCBD_corr) = NULL

# COntrib de chaque lakes à la BD_Total
SSi = rowSums(S)
df_LCBD_corr = as.data.frame(SSi / SS_Total)
colnames(df_LCBD_corr) = "LCBD_corr"
df_LCBD_corr$lake_id = as.factor(rownames(df_LCBD_corr))
rownames(df_LCBD_corr) = NULL


#####################################################################################################
# Export CSV ---------------------------------------------------------------------------------------
#####################################################################################################

data_python = env[,-c(38:46)] %>%
  left_join(summary_data_phyto[,-c(2:4,6)],by="lake_id") %>%
  left_join(summary_data_phyto_corr[,-2],by="lake_id") %>%
  left_join(summary_data_zoo[,-c(6:10)],by="lake_id") %>%
  
  left_join(df_LCBD[,c("lake_id","LCBD")],by="lake_id") %>%
  left_join(df_LCBD_corr[,c("lake_id","LCBD_corr")],by="lake_id") %>%
  
  mutate(log_LCBD = log(LCBD)) %>%
  mutate(log_LCBD_corr = log(LCBD_corr)) %>%
  dplyr::select(-c(Survey,secchi_bottom,Stratification,LCBD,LCBD_corr))

data_python = data_python[!is.na(data_python$shannon), ]
data_python = data_python[!is.na(data_python$eveness_piel), ]

data_python = data_python[!is.na(data_python$shannon_corr), ]
data_python = data_python[!is.na(data_python$eveness_piel_corr), ]

write.csv(data_python,"/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python/df_LP_NLA.csv",row.names=F)

#####################################################################################################
# Analyses ---------------------------------------------------------------------------------------
#####################################################################################################
# 
# # --- Div phyto ~ Predictors ---------------------------------------------------------------------------
# numeric_vars = data_python %>%
#   dplyr::select(where(is.numeric))
# 
# for (var in names(numeric_vars)) {
#   gg = ggplot(data_python, aes(x = prev_Mixo, y = .data[[var]])) +
#     geom_point(alpha = 0.5) +
#     geom_smooth(method = "gam", se = T, color = "red") + # ajuste un gam sur chaque var
#     labs(x = "prev_Mixo",
#          y = var) +
#     theme_minimal()
# 
#   print(gg)
# }
# 
# # --- Div zoo ~ Predictors ---------------------------------------------------------------------------
# data_zoo_test = env[,-c(38:46)] %>%
#   left_join(summary_data_zoo,by="lake_id") %>%
#   dplyr::select(-c(Survey,secchi_bottom,Stratification))
# 
# 
# numeric_vars = data_zoo_test %>%
#   dplyr::select(where(is.numeric))
# 
# for (var in names(numeric_vars)) {
#   gg = ggplot(numeric_vars, aes(x = prev_Mixo, y = .data[[var]])) +
#     geom_point(alpha = 0.5) +
#     geom_smooth(method = "gam", se = T, color = "red") + # ajuste un gam sur chaque var
#     labs(x = "prev_Mixo",
#          y = var) +
#     theme_minimal()
# 
#   print(gg)
# }






#####################################################################################################
# Autres relations verif ----------------------------------------------------------------------------
#####################################################################################################


# -----------------------------------------------------------------------------------------------
# Relation prev_Cyano ~ Diversité TOTALE                                                                             
# -----------------------------------------------------------------------------------------------
summary_data_cyano = data_phyto %>%
  rename("genus" = "Target_taxon",
         "biovol_tot" = "total.biov") %>%
  mutate(strat = ifelse(genus %in% all_mixo, "Mixotroph", "Autotroph")) %>%
  filter(biovol > 0) %>%
  group_by(lake_id)%>%
  summarize(
    Survey = first(Survey),
    biovol_tot_um3.mL = first(biovol_tot),
    biovol_mixo_um3.mL = sum(biovol[strat == "Mixotroph"]),
    biovol_cyano_um3.mL = sum(biovol[classic.group == "CYANOBACTERIA"]),
    rich_genus = n_distinct(genus),
    rich_genus_mixo = sum(ifelse(strat=="Mixotroph",1,0)),
    shannon = -sum((biovol/biovol_tot_um3.mL) * log(biovol/biovol_tot_um3.mL) ),
    op_simpson = 1-sum((biovol/biovol_tot_um3.mL)**2),
  ) %>%
  ungroup() %>%
  mutate(eveness_piel = shannon/log(rich_genus),
         prev_Mixo = biovol_mixo_um3.mL/biovol_tot_um3.mL*100,
         prev_Cyano = biovol_cyano_um3.mL/biovol_tot_um3.mL*100) 


plot(x=summary_data_cyano$prev_Cyano, y= summary_data_cyano$rich_genus)
plot(x=summary_data_cyano$prev_Cyano, y= summary_data_cyano$shannon)
plot(x=summary_data_cyano$prev_Cyano, y= summary_data_cyano$op_simpson)
plot(x=summary_data_cyano$prev_Cyano, y= summary_data_cyano$eveness_piel)

numeric_vars = summary_data_cyano %>%
  dplyr::select(where(is.numeric))

for (var in names(numeric_vars)) {
  gg = ggplot(numeric_vars, aes(x = prev_Cyano, y = .data[[var]])) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "gam", se = T, color = "red") + # ajuste un gam sur chaque var
    labs(x = "prev_cyano",
         y = var) +
    theme_minimal()

  print(gg)
}

# -----------------------------------------------------------------------------------------------
# relation diversite prev_Mixo ~ diversite Autotrophs sans cyanobacteries                                                                             
# -----------------------------------------------------------------------------------------------
summary_data_phyto_corr2 = data_phyto %>%
  rename("genus" = "Target_taxon",
         "biovol_tot" = "total.biov") %>%
  dplyr::select(-c("biovol_tot"))  %>%
  mutate(strat = ifelse(genus %in% all_mixo, "Mixotroph", "Autotroph")) %>%
  filter(biovol > 0) %>%
  filter(strat=="Autotroph",
         classic.group!="CYANOBACTERIA") %>%
  group_by(lake_id) %>%
  summarize(
    biovol_tot_um3.mL = sum(biovol),
    
    rich_genus_corr = n_distinct(genus),
    shannon_corr = -sum((biovol/biovol_tot_um3.mL) * log(biovol/biovol_tot_um3.mL) ),
    op_simpson_corr = 1-sum((biovol/biovol_tot_um3.mL)**2),
  ) %>%
  ungroup() %>%
  mutate(eveness_piel_corr = shannon_corr/log(rich_genus_corr)) %>% 
  left_join(summary_data_phyto[,c(1,10)],by="lake_id")

plot(x=summary_data_phyto_corr2$prev_Mixo, y= summary_data_phyto_corr2$rich_genus_corr)
plot(x=summary_data_phyto_corr2$prev_Mixo, y= summary_data_phyto_corr2$shannon_corr)
plot(x=summary_data_phyto_corr2$prev_Mixo, y= summary_data_phyto_corr2$op_simpson_corr)
plot(x=summary_data_phyto_corr2$prev_Mixo, y= summary_data_phyto_corr2$eveness_piel_corr)

numeric_vars1 = summary_data_phyto_corr %>%
  left_join(summary_data_phyto[,c(1,10)],by="lake_id") %>%
  dplyr::select(where(is.numeric),
                -c(biovol_tot_um3.mL)) 

numeric_vars2 = summary_data_phyto_corr2 %>%
  dplyr::select(where(is.numeric),
                -c(biovol_tot_um3.mL)) 

for (var in names(numeric_vars1)) {
  if (var != "prev_Mixo"){
    
    if(var == "rich_genus_corr"){method_mod = "glm"}else{method_mod="gam"}
    
    gg1 = ggplot(numeric_vars1, aes(x = prev_Mixo, y = .data[[var]])) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = method_mod, se = T, color = "red") + # ajuste un gam sur chaque var
      labs(x = "prev_Mixo",
           y = paste("Div ", var," Autotrophes")) +
      theme_minimal()
  
    gg2 = ggplot(numeric_vars2, aes(x = prev_Mixo, y = .data[[var]])) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = method_mod, se = T, color = "red") + # ajuste un gam sur chaque var
      labs(x = "prev_Mixo",
           y = paste("Div ", var," Autotrophes SANS CYANO")) +
      theme_minimal()
    
    grid.arrange(gg1,gg2)
  }
}





#####################################################################################################
# DATA ZOO MANIPS -----------------------------------------------------------------------------------
#####################################################################################################

data_zoo = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/2_LP_NLA/ZooLakePulseRaw_ALLfinal_grouping2017_2018_2019.csv",sep=";") %>%
  rename("lake_id"="Lake_ID")
# rownames(data_zoo) = data_zoo$lake_id

H = diversity(as.matrix(data_zoo[,-1]), index="shannon")
one_minus_D = diversity(as.matrix(data_zoo[,-1]), index="simpson") # return 1-D

summary_data_zoo_spe = data_zoo %>%
  rowwise() %>%
  mutate(
    biovol_tot_zoo_um3.mL = sum(c_across(-lake_id), na.rm = TRUE),
    rich_spe = sum(c_across(!any_of("lake_id")) > 0)
  ) %>%
  ungroup() %>%
  mutate(shannon = H,
         op_simpson = one_minus_D, 
         eveness_piel = H/log(rich_spe)) %>%
  left_join(summary_data_phyto[,c(1,10)],by="lake_id") %>%
  dplyr::select(lake_id, biovol_tot_zoo_um3.mL,
                rich_spe, shannon, op_simpson, eveness_piel, prev_Mixo)


plot(x=summary_data_zoo_spe$prev_Mixo, y= summary_data_zoo_spe$rich_spe)
plot(x=summary_data_zoo_spe$prev_Mixo, y= summary_data_zoo_spe$shannon)
plot(x=summary_data_zoo_spe$prev_Mixo, y= summary_data_zoo_spe$op_simpson)
plot(x=summary_data_zoo_spe$prev_Mixo, y= summary_data_zoo_spe$eveness_piel)

numeric_vars = summary_data_zoo_spe %>%
  dplyr::select(where(is.numeric))

for (var in names(numeric_vars)) {
  gg = ggplot(numeric_vars, aes(x = prev_Mixo, y = .data[[var]])) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "gam", se = T, color = "red") + # ajuste un gam sur chaque var
    labs(x = "prev_Mixo",
         y = paste(var," zoo")) +
    theme_minimal()
  
  print(gg)
}





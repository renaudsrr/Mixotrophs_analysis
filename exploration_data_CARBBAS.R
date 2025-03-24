# Script by SERRE Renaud
# Internship UQAM-GRIL under supervision of BEISNER Beatrix

library(dplyr)
library(openxlsx)
library(tidyverse)
library(gridExtra)
library(corrplot)

rm(list=ls()) 
graphics.off() 
cat("\014")

setwd("/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/Data")

###########################################################################
# Data --------------------------------------------------------------------
###########################################################################

# phyto (genus)  ---------------------------------------------------------------------------

data_phyto_raw = read.csv("CARBBAS/data_phyto_genus.csv",sep=";",h=T,dec=",")
data_phyto = read.csv("CARBBAS/data_phyto_genus.csv",sep=";",h=F,dec=",") # creation data_phyto pour plus de lisibilite
data_phyto = data_phyto[-c(1,3),-1]
data_phyto[1,1] = "Sample"
colnames(data_phyto) = as.character(unlist(data_phyto[1,]))
data_phyto = data_phyto[-1,] 

data_phyto[, -1] = lapply(data_phyto[, -1], function(x) { # forcer conversion "," en ".", car pas fait par read.csv !
  if (is.character(x)) as.numeric(gsub(",", ".", x)) else x # code chatGPT
})

# calcul somme biomass
biomass_tot = rowSums(data_phyto[,-1]) 
summary_data_phyto = data_phyto %>%
  mutate("biomass_tot" = biomass_tot) %>% 
  select(Sample, biomass_tot)

# calcul richesse spé
summary_data_phyto = summary_data_phyto %>%
  mutate(rich_genus = rowSums(ifelse(data_phyto[,-1]!= 0,1,0)))

# zoo (sps/genus) -----------------------------------------------------------------------------

data_zoo  = read.csv('CARBBAS/data_zoo_sp.csv',sep=";",h=T, dec=",") 
data_zoo = data_zoo %>% 
  select(-c("X","X.1","X.2")) %>%
  rename(Sample = ALL.LAKES..)

# calcul somme biomass
biomass_tot = rowSums(data_zoo[,-1]) 
summary_data_zoo = data_zoo %>% 
  mutate("biomass_tot" = biomass_tot) %>% 
  select(Sample, biomass_tot) 

# calcul richesse spé
summary_data_zoo = summary_data_zoo %>%
  mutate(rich_spe = rowSums(ifelse(data_zoo[,-1]!= 0,1,0))) %>% 
  select(Sample, biomass_tot,rich_spe, everything()) 

# env --------------------------------------------------------------------------------------------

env = list("carbon_metabol" = read.xlsx("CARBBAS/EnvData/EnvData.xlsx", sheet=1),
           "biol" = read.xlsx("CARBBAS/EnvData/EnvData.xlsx", sheet=2),
           "phy" = read.xlsx("CARBBAS/EnvData/EnvData.xlsx", sheet=3),
           "chim" = read.xlsx("CARBBAS/EnvData/EnvData.xlsx", sheet=4),
           "GPP" = read.xlsx("CARBBAS/EnvData/EnvData.xlsx", sheet=5))
colnames(env$carbon_metabol)[6] = "Sample"
colnames(env$phy)[5] = "Sample"
colnames(env$chim)[5] = "Sample"

data_python  = env$carbon_metabol %>%
  select(-c("region","lake","yr","doy")) %>%
  right_join(env$biol,by=c("Sample")) %>% select(-c("region.x")) %>% rename("region"="region.y") %>%
  select("Sample","region","lake","yr","doy",everything()) 
data_python = env$chim %>%
  select(-c("region","lake","yr","doy")) %>%
  left_join(data_python,by=c("Sample")) %>%
  select("Sample","region","lake","yr","doy",everything()) 
data_python = env$phy %>%
  select(-c("region","lake","yr","doy")) %>%
  left_join(data_python,by=c("Sample")) %>%
  select("Sample","region","lake","yr","doy",everything()) 
data_python = data_python %>%
  mutate("zmix/zmax"=zmix/zmax) %>%
  mutate("type_prof"=ifelse(`zmix/zmax`<0.9,"stratifie","polymictique")) %>%
  select(-c("pctdo","pctdo.cor")) %>% 
  full_join(summary_data_zoo[,c(1:2)],by=c("Sample")) %>% 
  rename("biomass_tot_zoo"="biomass_tot",
        "r.ugcld"  = "r.ugcld.(bottle.R)")
data_python = data_python[c(1:69),]

# setdiff(summary_data_phyto$Sample,data_python$Sample)

 # Categorisation genus => Creation dataset info_genus_zoo for summary of all informations ---------

names_phyto = read.csv("Mixotroph_strat/Names_phyto.csv",sep=";",h=T)
plancton_strat = read.csv("Mixotroph_strat/NanoplanktonNutritionStrategies.csv",sep=";",h=T) %>% select(c(1:9),13)

old_strat_genus = data.frame(t(data_phyto_raw[1,-c(1,2)]),
                      t(data_phyto_raw[2,-c(1,2)])) 
colnames(old_strat_genus) = c("Abbreviation", "strategy_old")

info_genus_zoo = names_phyto %>% left_join(old_strat_genus, by = "Abbreviation")
info_genus_zoo = info_genus_zoo %>% mutate(Genus = substr(Genus, 1, nchar(Genus) - 4))
info_genus_zoo = info_genus_zoo %>% right_join(plancton_strat, by="Genus") %>%
  select(Genus,Abbreviation,strategy_old,Final_Nutrition_Strategy,everything()) %>% 
  mutate(strategy_old = ifelse(strategy_old == "Y","Mixotroph","Autotroph"))

###########################################################################
# Metrics richness --------------------------------------------------------------------
###########################################################################
source("/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/fct_div.R")

# + ajout nombre mixotrophe dans data_phyto2
all_mixo = info_genus_zoo$Abbreviation[!is.na(info_genus_zoo$Abbreviation) & info_genus_zoo$Final_Nutrition_Strategy == "Mixotroph"]
length(all_mixo)

calc_infos_mixo = function(data){
  nb_mixo = numeric(length(data$Sample))
  biomass_tot_mixo = numeric(length(data$Sample))
  Genus = colnames(data)[-1]  # all col genus
  
  for (lake in 1:length(data$Sample)){
    sum_mixo = 0
    biomass_mixo = 0
    
    for (genus in Genus){
      biomass_i = data[lake, genus]  
      
      if (biomass_i != 0 & genus %in% all_mixo){  
        sum_mixo = sum_mixo + 1
        biomass_mixo = biomass_mixo + biomass_i
      }
    }
    nb_mixo[lake] = sum_mixo
    biomass_tot_mixo[lake] = biomass_mixo
  }
  infos_mixo  = list(nb_mixo,biomass_tot_mixo)
  return(infos_mixo)
}

summary_data_phyto = summary_data_phyto %>%
  mutate(H = Shannon(data_phyto)) %>%
  mutate("1-D" = Simpson(data_phyto)) %>%
  mutate(J = H/log(rich_genus)) %>% # ajout equitabilite Pielou
  mutate(nb_genus_mixo = calc_infos_mixo(data_phyto)[[1]]) %>% # nb mixotrophe
  mutate(biomass_mixo = calc_infos_mixo(data_phyto)[[2]]) %>% # biomass mixotrophe
  mutate(prev_Mixo = biomass_mixo/biomass_tot*100) %>%
  select(Sample, biomass_tot,biomass_mixo,rich_genus, nb_genus_mixo, H,"1-D", J, prev_Mixo)

index_mixo = which(colnames(data_phyto) %in% all_mixo) ; print(index_mixo)
colnames(data_phyto[index_mixo]) # petite verif

###########################################################################
# Analyse -----------------------------------------------------------------
###########################################################################

# phyto (genus) : data_phyto -----------------------------------------------------------------

gg1 = ggplot(data=summary_data_phyto,aes(x=Sample))+
  geom_col(aes(y=rich_genus),alpha=0.8)+
  theme(axis.text.x = element_text(angle = 90))
gg2 = ggplot(data=summary_data_phyto,aes(x=Sample))+
  geom_col(aes(y=biomass_tot),alpha=0.8,fill="grey")+
  theme(axis.text.x = element_text(angle = 90))
gg3 = ggplot(summary_data_phyto,aes(x=Sample))+
  geom_col(aes(y=prev_Mixo),alpha=0.7)
grid.arrange(gg1,gg2,gg3,ncol=2,nrow=2)

# zoo (sps/genus) : data_zoo -----------------------------------------------------------------

gg1 = ggplot(data=summary_data_zoo,aes(x=Sample))+
  geom_col(aes(y=rich_spe),alpha=0.8)+
  theme(axis.text.x = element_text(angle = 90))
gg2 = ggplot(data=summary_data_zoo,aes(x=Sample))+
  geom_col(aes(y=biomass_tot),alpha=0.8,fill="grey")+
  theme(axis.text.x = element_text(angle = 90))

grid.arrange(gg1,gg2,ncol=2)

# env : env -----------------------------------------------------------------

corrplot(cor(env$carbon_metabol[,7:10]))
corrplot(cor(env$carbon_metabol[,8:10]))
corrplot(cor(env$phy[,6:12]))
corrplot(cor(env$GPP[,c(4:20,22,24:27)]))

# Categorisation genus : info_genus_zoo -----------------------------------------------------------------

# ou infos_genus_zoo tout court pour graph sur l'ensemble des genus
ggplot(data = info_genus_zoo[1:74,], aes(x=Classic.group.name, fill=Final_Nutrition_Strategy)) +
  geom_bar() +
  labs(x="Taxon", y="nombre d'espèces dans jeu de donnée", fill="Strategy")+
  theme(axis.text.x = element_text(angle = 45,hjust=1))

###########################################################################
# export CSV --------------------------------------------------------------
###########################################################################

data_python = data_python %>%
  left_join(summary_data_phyto[,c(1,4,6:9)],by="Sample")
data_python = data_python[match(summary_data_phyto$Sample, data_python$Sample), ]

write.csv(data_phyto, "new_csv/data_count/data_phyto_CARBBAS.csv", row.names=F)
write.csv(data_zoo, "new_csv/data_count/data_zoo_CARBBAS.csv", row.names=F)
write.csv(summary_data_phyto, "new_csv/summary_data/summary_data_phyto_CARBBAS.csv", row.names=F)
write.csv(summary_data_zoo, "new_csv/summary_data/summary_data_zoo_CARBBAS.csv", row.names=F)
write.csv(info_genus_zoo, "new_csv/info_genus.csv", row.names=F)
write.csv(data_python,"new_csv/data_python/data_python_CARBBAS.csv",row.names=F)

###########################################################################
# verification --------------------------------------------------------------
###########################################################################
# relation entre preMixo et diversité

plot(x = summary_data_phyto$prev_Mixo, y = summary_data_phyto$rich_genus)
plot(x = data_python$prev_Mixo, y = summary_data_phyto$rich_genus,col="red")

plot(x = summary_data_phyto$prev_Mixo, y = summary_data_phyto$H)
plot(x = summary_data_phyto$prev_Mixo, y = summary_data_phyto$`1-D`)
plot(x = summary_data_phyto$prev_Mixo, y = summary_data_phyto$J)
# tendance logistique puis décroissante ?



###########################################################################
# QQ infos supp -----------------------------------------------------------
###########################################################################

lakes_phyto = sort(unique(summary_data_phyto$Sample))
lakes_zoo = sort(unique(summary_data_zoo$Sample))
lakes_env = unique(env$carbon_metabol$Sample) # 1 à 5 lakes env indentiques

all_values = sort(unique(c(lakes_phyto, lakes_zoo, lakes_env)))

aligned_lakes_phyto = ifelse(all_values %in% lakes_phyto, all_values, NA)
aligned_lakes_zoo = ifelse(all_values %in% lakes_zoo, all_values, NA)
aligned_lakes_env = ifelse(all_values %in% lakes_env, all_values, NA)

df_diff_lakes = data.frame(lakes_env = aligned_lakes_env,
                           lakes_phyto = aligned_lakes_phyto, 
                           lakes_zoo = aligned_lakes_zoo)

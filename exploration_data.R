# Script by SERRE Renaud
# Internship UQAM-GRIL under supervision of BEISNER Beatrix

library(dplyr)
library(openxlsx)
library(tidyverse)
library(gridExtra)


rm(list=ls()) 
graphics.off() 
cat("\014")

setwd("/Users/renaudsrr/Desktop/STAGE_MTL/Data")

###########################################################################
# Data --------------------------------------------------------------------
###########################################################################

# phyto (genus)  ---------------------------------------------------------------------------

data_phyto_raw = read.csv("data_phyto_genus.csv",sep=";",h=T,dec=",")
data_phyto = read.csv("data_phyto_genus.csv",sep=";",h=F,dec=",") # creation data_phyto pour plus de lisibilite
data_phyto = data_phyto[-c(1,3),-1]
data_phyto[1,1] = "Lakes"
colnames(data_phyto) = as.character(unlist(data_phyto[1,]))
data_phyto = data_phyto[-1,] 

data_phyto[, -1] = lapply(data_phyto[, -1], function(x) { # forcer conversion "," en ".", car pas fait par read.csv !
  if (is.character(x)) as.numeric(gsub(",", ".", x)) else x # code chatGPT
})

# calcul somme biomass
biomass_tot = rowSums(data_phyto[,-1]) 
data_phyto = data_phyto %>%
  mutate("biomass_tot" = biomass_tot) %>% 
  select(Lakes, biomass_tot,everything())

# calcul richesse spé
data_phyto = data_phyto %>%
  mutate(rich_spe = rowSums(ifelse(.[,-c(1,2)]!= 0,1,0))) %>% 
  select(Lakes, biomass_tot,rich_spe, everything()) 


# zoo (sps/genus) -----------------------------------------------------------------------------

data_zoo  = read.csv('Data_zoo_sp.csv',sep=";",h=T, dec=",") 
data_zoo = data_zoo %>% 
  select(-c("X","X.1","X.2")) %>%
  rename(Lakes = ALL.LAKES..)

# calcul somme biomass
biomass_tot = rowSums(data_zoo[,-1]) 
data_zoo = data_zoo %>% 
  mutate("biomass_tot" = biomass_tot) %>% 
  select(Lakes, biomass_tot,everything()) 

# calcul richesse spé
data_zoo = data_zoo %>%
  mutate(rich_spe = rowSums(ifelse(.[,-c(1,2)]!= 0,1,0))) %>% 
  select(Lakes, biomass_tot,rich_spe, everything()) 

# env --------------------------------------------------------------------------------------------

env = list("carbon_metabol" = read.xlsx("EnvData.xlsx", sheet=1),
           "biol" = read.xlsx("EnvData.xlsx", sheet=2),
           "phy" = read.xlsx("EnvData.xlsx", sheet=3),
           "chim" = read.xlsx("EnvData.xlsx", sheet=4),
           "GPP" = read.xlsx("EnvData.xlsx", sheet=5))

# Categorisation genus => Creation dataset info_genus_zoo for summary of all informations ---------

names_phyto = read.csv("Names_phyto.csv",sep=";",h=T)
plancton_strat = read.csv("NanoplanktonNutritionStrategies.csv",sep=";",h=T) %>% select(c(1:9),13)

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

Shannon = function(data){
  H = numeric(length(data$Lakes)) 
  
  for (lake in 1:length(data$Lakes)){
    total_biomass = data$biomass_tot[lake]
    sum_H = 0  
    
    for (genus in 4:ncol(data)){  # Boucle genus
      biomass_i = data[lake, genus]  # Biomasse sp i
      
      if (biomass_i > 0){  # Eviter log(0)
        p_i = biomass_i / total_biomass
        sum_H = sum_H + (p_i * log(p_i))}
    }
    H[lake] = -sum_H
  }
  return(H)
}

Simpson = function(data){
  D = numeric(length(data$Lakes)) 
  
  for (lake in 1:length(data$Lakes)){
    total_biomass = data$biomass_tot[lake]
    sum_D = 0
    
    for (genus in 4:ncol(data)){  # Boucle genus
      biomass_i = data[lake, genus]  # Biomasse sp i
      
      if (biomass_i > 0){  # evite calcul + 0
        p_i = biomass_i / total_biomass
        sum_D = sum_D + p_i**2}
    }
    D[lake] = 1-sum_D # indide de Simpson complémentaire + intuitif
  }
  return(D)
}

data_phyto2 = data_phyto %>%
  mutate(H = Shannon(.)) %>%
  mutate("1-D" = Simpson(.)) %>%
  mutate(Equi = H/log(rich_spe)) %>%
  select(Lakes, biomass_tot,rich_spe, H,"1-D", Equi, everything())

data_zoo2 = data_zoo %>%
  mutate(H = Shannon(.)) %>%
  mutate("1-D" = Simpson(.)) %>%
  mutate(Equi = H/log(rich_spe)) %>%
  select(Lakes, biomass_tot,rich_spe, H,"1-D", Equi, everything())

###########################################################################
# Analyse -----------------------------------------------------------------
###########################################################################

# phyto (genus) : data_phyto -----------------------------------------------------------------

summary(data_phyto)

gg1 = ggplot(data=data_phyto,aes(x=Lakes))+
  geom_col(aes(y=rich_spe),alpha=0.8)+
  theme(axis.text.x = element_text(angle = 90))
gg2 = ggplot(data=data_phyto,aes(x=Lakes))+
  geom_col(aes(y=biomass_tot),alpha=0.8,fill="grey")+
  theme(axis.text.x = element_text(angle = 90))

grid.arrange(gg1,gg2,ncol=2)

# zoo (sps/genus) : data_zoo -----------------------------------------------------------------

gg1 = ggplot(data=data_zoo,aes(x=Lakes))+
  geom_col(aes(y=rich_spe),alpha=0.8)+
  theme(axis.text.x = element_text(angle = 90))
gg2 = ggplot(data=data_zoo,aes(x=Lakes))+
  geom_col(aes(y=biomass_tot),alpha=0.8,fill="grey")+
  theme(axis.text.x = element_text(angle = 90))

grid.arrange(gg1,gg2,ncol=2)

# env : env -----------------------------------------------------------------



# Categorisation genus : info_genus_zoo -----------------------------------------------------------------

# ou infos_genus_zoo tout court pour graph sur l'ensemble des genus
ggplot(data = info_genus_zoo[1:74,], aes(x=Classic.group.name, fill=Final_Nutrition_Strategy)) +
  geom_bar() +
  labs(x="Taxon", y="nombre d'espèces dans jeu de donnée", fill="Strategy")+
  theme(axis.text.x = element_text(angle = 45,hjust=1))




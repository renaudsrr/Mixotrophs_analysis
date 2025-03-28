# Script by SERRE Renaud
# Internship UQAM-GRIL under supervision of BEISNER Beatrix

library(dplyr)
library(openxlsx)
library(tidyverse)
library(gridExtra)


rm(list=ls()) 
graphics.off() 
cat("\014")

setwd("/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/Data/")

###########################################################################
# Import Data --------------------------------------------------------------------
###########################################################################

# sp/genus
data_phyto = read.csv("LP&NLA/HarmonizedPhyto_LakePulseNLA2017_19092022.csv")
# data_phyto2 = read.csv("nla-2017-zooplankton-count-data.csv") # uniquement NLA
data_zoo = read.csv("LP&NLA/ZooLakePulseRaw_ALLfinal_grouping2017_2018_2019.csv",sep=";") %>%
  rename("lake_id"="Lake_ID")


# env
env = read.csv("LP&NLA/HarmonizedPredictors_LakePulseNLA2017_21062022.csv") %>%
  rename("lake_id"="Lake_ID") %>% select(Survey,everything())

# strat
info_genus = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/Data/new_csv/info_genus.csv")
# plancton_strat = read.csv("MatMix_AnnexeL_NLA2017LP_05012023.csv")
# plancton_data = read.xlsx("NLA2017LakePulse-zooplankton-taxa-list-data-06062022.xlsx")

###########################################################################
# Data manipulation --------------------------------------------------------------------
###########################################################################

all_mixo = info_genus$Genus[!is.na(info_genus$Genus) & info_genus$Final_Nutrition_Strategy == "Mixotroph"]
length(all_mixo)

summary_data_phyto = data_phyto %>%
  rename("genus" = "Target_taxon",
         "biovol_tot" = "total.biov") %>%
  select(-c("classic.group"))  %>%
  mutate(strat = ifelse(genus %in% all_mixo, "Mixotroph", "Autotroph")) %>%
  group_by(lake_id)%>%
  summarize(
    Survey = first(Survey),
    biovol_tot = first(biovol_tot),
    biovol_mixo = sum(biovol[strat == "Mixotroph"]),
    rich_genus = length(lake_id),
    nb_genus_mixo = sum(ifelse(strat=="Mixotroph",1,0)),
    H = -sum((biovol/biovol_tot) * log(biovol/biovol_tot) ),
    "1-D" = sum((biovol/biovol_tot)**2),
  ) %>% 
  mutate("J" = H/log(rich_genus),
         "prev_Mixo" = biovol_mixo/biovol_tot*100) %>%
  select("Survey",everything())

data_python = env %>%
  full_join(summary_data_phyto,by="lake_id") %>% arrange("lake_id") %>%
  rename("Survey" = "Survey.x") %>% select(-c("Survey.y")) %>%
  mutate(secchi_bottom = ifelse(secchi_bottom=="Yes",1,0))
a = dim(data_python)[1]
data_python = data_python[!is.na(data_python$H), ]
b = dim(data_python)[1]
print(a-b)
data_python = data_python[!is.na(data_python$J), ]
c = dim(data_python)[1]
print(a-c)

###########################################################################
# Traitement des NA --------------------------------------------------------------
###########################################################################
colSums(is.na(data_python))
dim(data_python)
###########################################################################
# export CSV --------------------------------------------------------------
###########################################################################
  
write.csv(summary_data_phyto, "new_csv/summary_data/summary_data_phyto_LP_NLA.csv", row.names=F)
write.csv(data_python,"new_csv/data_python/data_python_LP_NLA.csv",row.names=F)

###########################################################################
# verification --------------------------------------------------------------
###########################################################################
# relation entre preMixo et diversité

plot(x = summary_data_phyto$prev_Mixo, y = summary_data_phyto$rich_genus)
plot(x = data_python$prev_Mixo, y = data_python$rich_genus,col="red")

plot(x = summary_data_phyto$prev_Mixo, y = summary_data_phyto$H)
plot(x = summary_data_phyto$prev_Mixo, y = summary_data_phyto$`1-D`)
plot(x = summary_data_phyto$prev_Mixo, y = summary_data_phyto$J)
# tendance en U ou inv U ??!!

###########################################################################
# QQ infos supp -----------------------------------------------------------
###########################################################################

setdiff(env$lake_id,summary_data_phyto$lake_id)
# [1] "06-456"  "07-002"  "07-014"  "07-033"  "06-133"  "08-164"  "08-166"  "08-169"  "08-179"  "08-180"  "08-183"  "08-197" 
# [13] "08-201"  "08-211"  "09-450"  "12-634"  "12-647"  "17-065"  "17-099"  "17-108"  "17-114"  "17-122"  "18-592"  "18-594" 
# [25] "2010814" "2010815" "2010866" "2010869" "2010968" "2010970" "2011067" "2011108" "2011163" "2011175" "2011176" "2011181"
# [37] "2011228" "2011273" "2011309" "2011322" "2011325" "2011381" "2011409" "2011437" "2011495" "2011507" "2011508"

setdiff(summary_data_phyto$lake_id,env$lake_id)

















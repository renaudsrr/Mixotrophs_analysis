# Script by SERRE Renaud
# Internship UQAM-GRIL under supervision of BEISNER Beatrix

library(dplyr)
library(corrplot)
library(gridExtra)

rm(list = ls())
graphics.off()
cat("\014") 

# --- Classification de base - by le Noach ---------------------------------------------------------------------------
info_genus = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/Mixotroph_strat/info_genus.csv")
info_genus = info_genus %>%
  rename(strat = Final_Nutrition_Strategy) %>%
  select("Empire","Kingdom","Phylum","Class","Order","Family","Genus","strat")

summary(info_genus)
colSums(is.na(info_genus))

# --- Classification Scott IISD-ELA ---------------------------------------------------------------------------

infos_genus = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/Mixotroph_strat/info_genus.csv")
sp_code = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/4_IISD-ELA/csv_traites/dt5.csv")
strat_ELA = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/4_IISD-ELA/csv_traites/dt6.csv")

strat_ELA = strat_ELA %>% # !!! MEAN BY GENRE !!!!
  left_join(sp_code,by="species_code") %>%
  rename(Genus = genus)
strat_ELA_genus = strat_ELA %>%
  filter(!is.na(Genus), !is.na(trophic_type)) %>%
  group_by(Genus) %>%
  count(trophic_type) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup()

infos_genus_plus = infos_genus %>%
  left_join(strat_ELA_genus, by = "Genus") %>%
  rename(strategy_new = Final_Nutrition_Strategy,
         strategy_ELA = trophic_type) %>%
  mutate(strategy_ELA = case_when(
    strategy_ELA == "a" ~ "Autotroph",
    strategy_ELA == "m" ~ "Mixotroph",
    strategy_ELA == "h" ~ "Heterotroph",
    TRUE ~ strategy_ELA)) %>%
  select(Empire, Kingdom, Phylum, Class, Order, Family, Genus,
         strategy_old, strategy_new, strategy_ELA)

table(infos_genus_plus$strategy_new, infos_genus_plus$strategy_ELA, useNA = "ifany")

infos_genus_plus %>%
  filter(!is.na(strategy_new), !is.na(strategy_ELA), strategy_new != strategy_ELA) %>%
  select(Genus, strategy_new, strategy_ELA)

write.csv(info_genus,"/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/taxo_info.csv")



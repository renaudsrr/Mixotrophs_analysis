# Script by SERRE Renaud
# Internship UQAM-GRIL under supervision of BEISNER Beatrix

rm(list=ls()) 
graphics.off() 
cat("\014")

library(dplyr)
library(ggplot2)
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

data_phyto = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/1_CARBBAS/Data_phyto_genus.csv", sep=";", h=F, dec=",")[-3, -1]
data_phyto[1,1] = "lake_id"
colnames(data_phyto) = sub("\\..*", "", data_phyto[1,])
data_phyto = data_phyto[-c(1:2), ]
colnames(data_phyto) = make.names(colnames(data_phyto), unique = TRUE)
data_phyto = data_phyto %>%
  mutate(across(-1, ~ as.numeric(gsub(",", ".", .)))) %>%
  mutate(Pediastrum = Pediastrum + Pediastrum.1,
         Tetraedron = Tetraedron + Tetraedron.1) %>%
  select(-c(Pediastrum.1,Tetraedron.1))

summary(data_phyto)
colSums(is.na(data_phyto))

# -----------------------------------------------------------------------------------------------
# Zoo                                                                           
# -----------------------------------------------------------------------------------------------

data_zoo  = read.csv('/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/1_CARBBAS/Data_zoo_sp.csv',sep=";",h=T, dec=",") %>%
  select(-c("X","X.1","X.2")) %>%
  rename(lake_id = ALL.LAKES..) 
colnames(data_zoo) = sub("\\..*", "", colnames(data_zoo))
colnames(data_zoo) = make.names(colnames(data_zoo), unique = TRUE)
data_zoo = data_zoo %>%
  mutate(Acant = Acant + Acant.1,
         Cyclo = Cyclo + Cyclo.1 + Cyclo.2,
         Chydo = Chydo + Chydo.1,
         Cerio = Cerio + Cerio.1 + Cerio.2 + Cerio.3,
         Daphnia = Daphnia + Daphn + Daphn.1 + Daphn.2 + Daphn.3 + Daphn.4 + Daphn.5 + Daphn.6 + Daphn.7 + Daphn.8 + Diaph.1,
         Diacy = Diacy + Diacy.1,
         Eubo = Eubo + Eubo.1 + Eubo.2,
         Eucyc = Eucyc + Eucyc.1,
         Lepto = Lepto + Lepto.1 + Lepto.2 + Lepto.3,
         Mesoc = Mesoc + Mesoc.1,
         Microc = Microc + Microc.1,
         Neobos = Neobos + Neobos.1 + Neobos.2,
         Skist = Skist + Skist.1,
         Sino = Sino + Sino.1) %>%
  select(-c(Acant.1,Cyclo.1,Cyclo.2,Chydo.1,Cerio.1,Cerio.2,Cerio.3,Daphnia,Daphn.1,Daphn.2,Daphn.3,Daphn.4,Daphn.5,Daphn.6,
            Daphn.7, Daphn.8, Diacy.1, Eubo.1 , Eubo.2,Eucyc.1,Lepto.1 , Lepto.2 , Lepto.3,Mesoc.1, 
            Microc.1,Neobos.1 , Neobos.2, Skist.1,Sino.1))

summary(data_zoo)
colSums(is.na(data_zoo))

# -----------------------------------------------------------------------------------------------
# Env                                                                             
# -----------------------------------------------------------------------------------------------

env_sheets = list("carbon_metabol" = read.xlsx("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/1_CARBBAS/EnvData.xlsx", sheet=1),
           "biol" = read.xlsx("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/1_CARBBAS/EnvData.xlsx", sheet=2),
           "phy" = read.xlsx("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/1_CARBBAS/EnvData.xlsx", sheet=3),
           "chim" = read.xlsx("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/1_CARBBAS/EnvData.xlsx", sheet=4))
colnames(env_sheets$carbon_metabol)[6] = "lake_id"
colnames(env_sheets$biol)[5] = "lake_id"
colnames(env_sheets$phy)[5] = "lake_id"
colnames(env_sheets$chim)[5] = "lake_id"

clean_env = list(
  env_sheets$carbon_metabol %>% select(-c("region")), # car region en double
  env_sheets$biol  %>% select(-c("region", "lake", "yr", "doy")),
  env_sheets$chim %>% select(-c("region", "lake", "yr", "doy")),
  env_sheets$phy %>% select(-c("region", "lake", "yr", "doy")))
env = reduce(clean_env, left_join, by = "lake_id")

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

calc_infos_mixo = function(data){
  rich_genus_mixo = numeric(nrow(data))
  biovol_tot_mixo = numeric(nrow(data))
  Genus = colnames(data)[-1]
  
  Genus_in_data = intersect(Genus, all_mixo)  
  
  for (lake in seq_len(nrow(data))){
    mixo_vals = as.numeric(unlist(data[lake, Genus_in_data]))
    rich_genus_mixo[lake] = sum(mixo_vals > 0)
    biovol_tot_mixo[lake] = sum(mixo_vals)
    }
  return(list(rich_genus_mixo, biovol_tot_mixo))
}

H = diversity(as.matrix(data_phyto[,-1]), index="shannon")
one_minus_D = diversity(as.matrix(data_phyto[,-1]), index="simpson") # return 1-D
biovol_mixo_um3.mL = calc_infos_mixo(data_phyto)[[2]]
rich_genus_mixo = calc_infos_mixo(data_phyto)[[1]]

summary_data_phyto = data_phyto %>%
  rowwise() %>%
  summarize(biovol_tot_um3.mL = sum(c_across(!any_of("lake_id"))),
            rich_genus = sum(c_across(!any_of("lake_id")) > 0)) %>%
  ungroup() %>%
  mutate(biovol_mixo_um3.mL = biovol_mixo_um3.mL, 
         rich_genus_mixo = rich_genus_mixo,
         shannon = H,
         op_simpson = one_minus_D, 
         eveness_piel = H/log(rich_genus),
         prev_Mixo = biovol_mixo_um3.mL / biovol_tot_um3.mL * 100) %>%
  mutate(lake_id = data_phyto$lake_id) %>%
  select(lake_id, biovol_tot_um3.mL, biovol_mixo_um3.mL,
         rich_genus, rich_genus_mixo, shannon, op_simpson, eveness_piel, prev_Mixo)

# -----------------------------------------------------------------------------------------------
# Zoo                                                                             
# -----------------------------------------------------------------------------------------------

# Classification en 3 groupes : Cladocera, Copepoda et OTHER
cladocera = c("Alona", "Alonop", "Algla", "Bosmi", "Bosmina", "Chydo", "Cerio", "Daphn", "Eubos", "Eubo", "Eurycercus", "Holop", "Ilyocryptus", "Laton",
              "Lepto", "Macrothrix", "Moina", "Neobos", "Onycho", "Ophryo", "Ortho",
              "Polyp", "Sida", "Sinobos", "Sino", "Tropo", "Clado_imm")
copepoda = c("Acant", "Acrop", "Calanoida", "Cyclopoida", "Cyc", "Cyclo", "Diacyclops", "Diacy",
             "Eucyc", "Episc", "Harpacticoida", "Mesoc", "Mcycl", "Cal", "Paracy")
other = c("nauplii")

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
  left_join(summary_data_phyto[,c(1,9)],by="lake_id") %>%
  select(lake_id, biovol_tot_zoo_um3.mL, 
         biovol_cladocera_um3.mL, biovol_copepoda_um3.mL, biovol_other_um3.mL, 
         rich_genus, shannon, op_simpson, eveness_piel, prev_Mixo)


#####################################################################################################
# Export CSV ---------------------------------------------------------------------------------------
#####################################################################################################

data_python = env %>%
  left_join(summary_data_phyto[,c(1,4,6:9)],by="lake_id") %>%
  left_join(summary_data_zoo[,c(1:5)],by="lake_id") %>%
  select(-c(lake,region)) %>% select(lake_id,everything())

write.csv(data_python,"/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python/df_CARBBAS.csv",row.names=F)

#####################################################################################################
# Analyses ---------------------------------------------------------------------------------------
#####################################################################################################

# numeric_vars = data_python %>%
#   select(where(is.numeric))
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

###########################################################################
# QQ infos supp -----------------------------------------------------------
###########################################################################

# lakes_phyto = sort(unique(summary_data_phyto$Sample))
# lakes_zoo = sort(unique(summary_data_zoo$Sample))
# lakes_env = unique(env$carbon_metabol$Sample) # 1 à 5 lakes env indentiques
# 
# all_values = sort(unique(c(lakes_phyto, lakes_zoo, lakes_env)))
# 
# aligned_lakes_phyto = ifelse(all_values %in% lakes_phyto, all_values, NA)
# aligned_lakes_zoo = ifelse(all_values %in% lakes_zoo, all_values, NA)
# aligned_lakes_env = ifelse(all_values %in% lakes_env, all_values, NA)
# 
# df_diff_lakes = data.frame(lakes_env = aligned_lakes_env,
#                            lakes_phyto = aligned_lakes_phyto, 
#                            lakes_zoo = aligned_lakes_zoo)
#                            
# 


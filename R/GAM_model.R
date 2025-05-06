# Script by SERRE Renaud
# Internship UQAM-GRIL under supervision of BEISNER Beatrix

rm(list=ls()) 
graphics.off() 
cat("\014")

library(dplyr)
library(openxlsx)
library(gridExtra)
library(grid) 
library(ggplot2) # graph package
library(tinytex) # Pour la sortie pdf
library(corrplot)# Correlation matrix calculus
library(plot3D)# For 3D plot
library(MASS)# Perform NB GLM
library(DHARMa)# Model diagnosis
library(rcompanion)# Model pseudo R²
library(lattice)# multipanel graphics
library(mgcv)# GAM
library(effects)
library(caret)

# ----------- Choix de la variable réponse et du dataset -----------
instructions = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/CODE/Instructions/df_instructions.csv")
response_variable = instructions$choix_div
choix_dataset = instructions$choix_dataset
correction_div = instructions$correction_div

# --- DATA ---------------------------------------------------------------------------
setwd("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/sep_train_test")

X_train = read.csv("X_train.csv")
X_test = read.csv("X_test.csv")
y_train = read.csv("y_train.csv")
y_test = read.csv("y_test.csv")

colSums(is.na(X_train))

if (choix_dataset =="C"){
  X_train = X_train %>% dplyr::select(-c(yr))
  X_test = X_test %>% dplyr::select(-c(yr))
}

if (choix_dataset =="ELA"){
  vars_to_interp = c("color", "pH", "Chla_ug.L", "DO_mg.L", "DOC_mg.L", "TP_ug.L", "TN_ug.L",
                     "SO_mg.L", "MG_mg.L", "CL_mg.L", "K_mg.L", "COND_uS.cm")
  
  for (var in vars_to_interp) {
    X_train[[var]] = with(X_train, approx(doy, get(var), xout = doy)$y)
    X_test[[var]] = with(X_test, approx(doy, get(var), xout = doy)$y)
  }
  X_train = X_train %>% dplyr::select(-c(lat,long,area_m2,depth)) # trop peu de valeurs differentes pour le spline
  X_test = X_test %>% dplyr::select(-c(lat,long,area_m2,depth))
}

if (choix_dataset =="LN"){
  X_train = X_train %>% dplyr::select(-c(biovol_tot_zoo_um3.mL,biovol_cladocera_um3.mL,biovol_copepoda_um3.mL,biovol_other_um3.mL))
  X_test = X_test %>% dplyr::select(-c(biovol_tot_zoo_um3.mL,biovol_cladocera_um3.mL,biovol_copepoda_um3.mL,biovol_other_um3.mL))
}

mask_train = complete.cases(X_train)
mask_test  = complete.cases(X_test)
X_train = X_train[mask_train, ]
y_train = y_train[mask_train, , drop = FALSE]
X_test  = X_test[mask_test, ]
y_test  = y_test[mask_test, , drop = FALSE]

X_train = na.omit(X_train)
X_test = na.omit(X_test)
train_set = cbind(X_train, y_train)
test_set = cbind(X_test, y_test)

##################################################################################
# GAM model ---------------------------------------------------------------
##################################################################################

# --- provisoire ---------------------------------------------------------------------------
if (choix_dataset =="ELA"){
  train_set = train_set %>% dplyr::select(-c(biovol_tot_zoo_um3.mL,biovol_cladocera_um3.mL,biovol_copepoda_um3.mL,biovol_other_um3.mL,sech_m,temp_mean))
  X_train = X_train %>% dplyr::select(-c(biovol_tot_zoo_um3.mL,biovol_cladocera_um3.mL,biovol_copepoda_um3.mL,biovol_other_um3.mL,sech_m,temp_mean))
  }

pred_names = names(X_train)[names(X_train) != response_variable]
smoothers = paste0("s(", pred_names, ")")
formula_gam = as.formula(paste(response_variable, "~", paste(smoothers, collapse = " + ")))
# Attention ici pas étude des effets combines

# scat(link="identity") 
# Gamma(link = "log")
# tw(link="log")

t0 = Sys.time()
mod1 = gam(formula_gam,
           family = gaussian() , # pas gaussian car Y asymétrique mais marche aussi avec gaussien donc bizzare...
           data=train_set,
           method = "REML",
           select = T)
summary(mod1)
t1 = Sys.time()
print(t1 - t0)

# Model avec parcimonie ---------------------------------------------------

# formula_reduced = as.formula(
#   paste(response_variable, "~",
#         "s(lat) + s(long) + s(area) + s(zmax) +",
#         "s(Temperature_up) + s(pH_up) + s(TP) + s(color) +",
#         "s(Chl) + s(biomass_tot_zoo) + s(prev_Mixo)")
# )
# 
# mod1 = gam(formula_reduced, data = X_train)
# summary(mod1)

# Verification model ------------------------------------------------------

# plot(mod1)
plot(mod1,pages=1,residuals=TRUE)  ## show partial residuals
plot(mod1,pages=1,seWithMean=TRUE) ## `with intercept' CIs

# See 3D plot --------------------------------------------------------------
# 
# # prev_Mixo et long
# formula_3d = as.formula(paste(response_variable, "~ s(prev_Mixo) + s(long)"))
# M1b = gam(formula_3d, family = gaussian, data = train_set)
# par(mar = c(2, 2, 2, 2))
# vis.gam(M1b, theta = -25, color = "heat")
# 
# # prev_Mixo et biomass_tot_zoo
# formula_3d = as.formula(paste(response_variable, "~ s(prev_Mixo) + s(biomass_tot_zoo)"))
# M1b = gam(formula_3d, family = gaussian, data = train_set)
# par(mar = c(2, 2, 2, 2))
# vis.gam(M1b, theta = -25, color = "heat")


# VALIDATION MODELISATION -------------------------------------------------

# éhchantillon independant ? => sinon GAMM
# ajouter interactions entre covariables ?

# >0.9 => a éliminer
# concurvity(mod1, full = TRUE) # verif Absence de concurvité forte (colinéarité non linéaire) 
# 
# gam.check(mod1) # vérif résidus aléatoires (homosc et normalité)
# 
# influ_values = influence.gam(mod1)  # ou plot(gam.influence(mod)) # Verification individus influents 
# plot(influ_values)

# Prediction --------------------------------------------------------------
preds = predict(mod1, newdata = X_test, se.fit = TRUE,type="response")
# preds = predict(mod1, newdata = read.csv("/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/Data/new_csv/data_python/data_python_CARBBAS_harmz.csv"), 
#                 se.fit = TRUE)

test_set$fit = preds$fit

ggplot(test_set, aes_string(x = "prev_Mixo")) +
  geom_point(aes_string(y = response_variable), alpha = 0.5, color = "dodgerblue3") +
  geom_point(aes(y = fit), alpha = 0.5, color = "red") +
  theme_minimal() +
  labs(title = paste("Prédiction GAM pour", response_variable))
# coef(mod1)

res = summary(mod1)
r2_adj_model = res$r.sq

# Export CSV : prev_Mixo + fit
if ("prev_Mixo" %in% names(test_set)) {
  export_df = data.frame(
    prev_Mixo = test_set$prev_Mixo,
    fit = test_set$fit,
    r2_adj = r2_adj_model
  )
  if (correction_div=="N"){
    filename = paste0("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python/prediction_fit_GAM_all_div/fit_", 
                    response_variable, "_sur_",choix_dataset,"_avec_corr_",correction_div,".csv")}
  if (correction_div=="Y"){
    filename = paste0("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python/prediction_fit_GAM_corr_div/fit_", 
                      response_variable, "_sur_",choix_dataset,"_avec_corr_",correction_div,".csv")}
  
  write.csv(export_df, filename, row.names = FALSE)
}






















# -----------------------------------------------------------------------------------------------
# ANALYSE DATASET SI BESOIN                                                                             
# -----------------------------------------------------------------------------------------------

# Y : H
# Xs : le reste
# 
# # Distrib Y -----------------------------------------------------
# par(mfrow=c(2,2))
# boxplot(X_train[[response_variable]])
# dotchart(X_train[[response_variable]]) # cleveland plot
# hist(X_train[[response_variable]])
# qqnorm(X_train[[response_variable]])
# qqline(X_train[[response_variable]])
# # non normalité !
# 
# # Distrib cov Xs ----------------------------------------------------------
# Xs = X_train %>% dplyr::select(-!!sym(response_variable))
# # pairs(Xs)
# 
# # Relation Y~Xs -----------------------------------------------------------
# for (i in (1:length(Xs))){
#   X = Xs[,i]
#   name = names(Xs[i])
#   plot(X_train[[response_variable]]~X,pch=16,col='dodgerblue3',xlab=name,ylab='H')
# }
# 
# # Colinearite Xs ----------------------------------------------------------
# M = cor(Xs)
# par(mfrow = c(1, 1))
# corrplot.mixed(M,upper="square",lower.col="black",tl.col="black",cl.cex=0.8,
#               tl.cex=0.7, number.cex=0.8)
# 
# # which((M > 0.7 | M < -0.7) & M!=1)








"""
Script building a model [...]

Author : SERRE Renaud
Creation 12.03.2025

"""

# =============================================================================
# Packages
# =============================================================================

import os
import sys
import time
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from datetime import datetime
import subprocess

from scipy.stats import uniform,stats
from statsmodels.stats.outliers_influence import OLSInfluence

import xgboost as xgb
from xgboost import XGBRegressor

import shap
import optuna
optuna.logging.set_verbosity(optuna.logging.WARNING)

from sklearn.model_selection import train_test_split, GridSearchCV, RandomizedSearchCV, cross_val_score
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from sklearn.inspection import PartialDependenceDisplay
from sklearn.preprocessing import StandardScaler

from open_interface import open_interface

# =======================================================================================================================
# =======================================================================================================================
# PARAM SCRIPT 
# =======================================================================================================================
# =======================================================================================================================

params = open_interface()

choix_dataset = params["choix_dataset"]
choix_div = params["choix_div"]

# refaire_optimisation = params["refaire_optimisation"]
# input_op = params["input_op"]

test_size = params["test_size"]
harmonisation_dataset = params["harmonisation_dataset"]
run_gam_model = params["run_gam_model"]
correction_div = params["correction_div"]

learning_rate = params["learning_rate"]
max_depth = params["max_depth"]
n_estimators = params["n_estimators"]
colsample_bytree = params["colsample_bytree"]
subsample = params["subsample"]
reg_alpha = params["reg_alpha"]
reg_lambda = params["reg_lambda"]
min_child_weight = params["min_child_weight"]
gamma = params["gamma"]

if correction_div=="N":
    fig_dir = f"/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/figures/all_div/{choix_div}"
else:
    fig_dir = f"/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/figures/corr_div/{choix_div}"

# VERS SCRIPT R
df_instructions = pd.DataFrame([{"choix_dataset": choix_dataset, "choix_div": choix_div, "correction_div" : correction_div}])
df_instructions.to_csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/CODE/Instructions/df_instructions.csv")

# =============================================================================
# Script R 
# =============================================================================

script_CARBBAS = "/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/CODE/R/1_Processing_data/1_Processing_CARBBAS.R"
script_LP_NLA = "/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/CODE/R/1_Processing_data/2_Processing_LP_NLA.R"
script_IISD_ELA = "/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/CODE/R/1_Processing_data/3_Processing_IISD_ELA.R"
script_harmz = "/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/CODE/R/Harmonisation_data.R"
script_GAM = "/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/CODE/R/GAM_model.R"

# =============================================================================
# Data
# =============================================================================

if harmonisation_dataset == "N":
    if choix_dataset == "C":
        print("Run script CARBBAS...")
        subprocess.run(["Rscript", script_CARBBAS], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)   

        data_python = pd.read_csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python/df_CARBBAS.csv",index_col=0)
        VA_expl = data_python.drop(columns=["rich_genus","shannon","op_simpson","eveness_piel"])

    if choix_dataset == "LN":
        print("Run script LP_NLA...")
        subprocess.run(["Rscript", script_LP_NLA], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

        data_python = pd.read_csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python/df_LP_NLA.csv",index_col=0)
        VA_expl = data_python.drop(columns=["rich_genus","shannon","op_simpson","eveness_piel","log_LCBD","rich_genus_corr","shannon_corr","op_simpson_corr","eveness_piel_corr","log_LCBD_corr"])

    if choix_dataset =="ELA":
        print("Run script ELA...")
        subprocess.run(["Rscript", script_IISD_ELA], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)   

        data_python = pd.read_csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python/df_IIDS_ELA.csv",index_col=0)
        VA_expl = data_python.drop(columns=["rich_genus","shannon","op_simpson","eveness_piel","rich_genus_corr","shannon_corr","op_simpson_corr","eveness_piel_corr"])

else:
    print("Run script harmonisation...")
    subprocess.run(["Rscript", script_harmz], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    if choix_dataset == "C":
        print("Run script CARBBAS...")
        subprocess.run(["Rscript", script_CARBBAS], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)   

        data_python = pd.read_csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python/1_harmz_CARBBAS.csv",index_col=0)
        VA_expl = data_python.drop(columns=["rich_genus","shannon","op_simpson","eveness_piel"])

    if choix_dataset == "LN":
        print("Run script LP_NLA...")
        subprocess.run(["Rscript", script_LP_NLA], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

        data_python = pd.read_csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python/2_harmz_LP_NLA.csv",index_col=0)
        VA_expl = data_python.drop(columns=["rich_genus","shannon","op_simpson","eveness_piel","log_LCBD","rich_genus_corr","shannon_corr","op_simpson_corr","eveness_piel_corr","log_LCBD_corr"])

    if choix_dataset =="ELA":
        print("Run script ELA...")
        subprocess.run(["Rscript", script_IISD_ELA], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)   

        data_python = pd.read_csv("/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python/3_harmz_IISD_ELA.csv",index_col=0)
        VA_expl = data_python.drop(columns=["rich_genus","shannon","op_simpson","eveness_piel","rich_genus_corr","shannon_corr","op_simpson_corr","eveness_piel_corr"])

# =============================================================================
# Extraction data pour ML
# =============================================================================

if correction_div=="Y":
    y = data_python[f"{choix_div}_corr"] # Variable réponse
else:
    y = data_python[choix_div]  

X_train, X_test, y_train, y_test = train_test_split(VA_expl, y, test_size=test_size,random_state = 42) # Séparation avec jeu learning et test, ajout random_state = 42 si besoin pour cadrer résultats predi

print("var train : ", np.var(y_train), " // var test : ", np.var(y_test))

export_dir = "/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/sep_train_test"
X_train.to_csv(os.path.join(export_dir, "X_train.csv"), index=False)
X_test.to_csv(os.path.join(export_dir, "X_test.csv"), index=False)

# Si y_train et y_test sont des Series, on peut les convertir en DataFrame pour avoir des entêtes
y_train.to_frame(name=choix_div).to_csv(os.path.join(export_dir, "y_train.csv"), index=False)
y_test.to_frame(name=choix_div).to_csv(os.path.join(export_dir, "y_test.csv"), index=False)

# =============================================================================
# RUN GAM
# =============================================================================

if correction_div=="N": 
    file_fit_GAM = f"/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python/prediction_fit_GAM_all_div/fit_{choix_div}_sur_{choix_dataset}_avec_corr_{correction_div}.csv"
else:
    file_fit_GAM = f"/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/DATA/0_new_csv/data_python/prediction_fit_GAM_corr_div/fit_{choix_div}_sur_{choix_dataset}_avec_corr_{correction_div}.csv"

if run_gam_model == "Y":
    print("Run GAM model...")
    start = time.time()
    subprocess.run(["Rscript", script_GAM], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    end = time.time()
    print(f"Durée : {end - start:.4f} secondes")
    print("\n")

fit_GAM = pd.read_csv(file_fit_GAM)

# =============================================================================
# Exploration data
# =============================================================================

plt.figure(figsize=(10, 8))

plt.scatter(X_train["prev_Mixo"],y_train,color="black",alpha=0.7,label="train")
plt.scatter(X_test["prev_Mixo"],y_test,color="black",alpha=0.3,label="test")
plt.xlabel("Prévalence Mixotrophs")
plt.ylabel(f"Diversite : {choix_div}")

plt.legend()
plt.tight_layout()
fig_name = f"2_{choix_div}_en_fct_prev_Mixo"
plt.savefig(f"{fig_dir}/{fig_name}", dpi=300, bbox_inches="tight")  
plt.show()

# =============================================================================
# Optimisation hyperparamètres
# =============================================================================

# if refaire_optimisation == "Y":

#     # Bayesien Optuna ==================================================================
#     if int(input_op)>0:

#         print("opti Op en cours...")

#         def objective(trial):
#             params = {
#                 'n_estimators': trial.suggest_int('n_estimators', 50, 2000),
#                 'max_depth': trial.suggest_int('max_depth', 3, 10),
#                 'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.3, log=True),
#                 'subsample': trial.suggest_float('subsample', 0.5, 1.0),
#                 'colsample_bytree': trial.suggest_float('colsample_bytree', 0.5, 1.0),
#                 'reg_lambda': trial.suggest_float('reg_lambda', 0.0, 1.0),
#                 'reg_alpha': trial.suggest_float('reg_alpha', 0.0, 1.0),
#                 'gamma' : trial.suggest_float('gamma',0.0,0.5)
#             }
            
#             xgb_model_Op = XGBRegressor(**params)
            
#             score = cross_val_score(xgb_model_Op, X_train, y_train, cv=5, scoring="r2",n_jobs=-1)
#             return score.mean()

#         study = optuna.create_study(direction='maximize')
#         study.optimize(objective, n_trials=input_op)

#         hyperparams_Op = pd.DataFrame.from_dict(study.best_params, orient="index", columns=["Opti_Op"])
#         hyperparams_Rd = pd.DataFrame(np.zeros((len(study.best_params), 1)), columns=["Opti_Rd"])

#         if int(input_rd) > 0:
#             hyperparams_Rd = pd.DataFrame.from_dict(random_search.best_params_, orient="index", columns=["Opti_Rd"])
#             all_optimisation = pd.merge(hyperparams_Rd, hyperparams_Op, left_index=True, right_index=True, how="outer")
#         else:
#             all_optimisation = hyperparams_Op

#         ordre_hyperparams = ["learning_rate", "max_depth", "n_estimators", "colsample_bytree", "subsample", "reg_alpha", "reg_lambda", "gamma"]
#         all_optimisation = all_optimisation.reindex(ordre_hyperparams)

#         all_optimisation.to_csv(f"/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/Data/new_csv/optimisation_hyperparms/all_optimisation_{choix_div}.csv")

#         print("\n")
#         print(f"Optimisation hyperparams pour y = {choix_div}")
#         print(all_optimisation)
#         print(f"Meilleur R² Optuna : {study.best_value:.4f}")
#         print("\n")

# =============================================================================
# Creation model
# =============================================================================

xgb_model_1 = XGBRegressor(
    random_state=42,
    learning_rate=learning_rate,
    max_depth=max_depth,
    n_estimators=n_estimators,
    colsample_bytree=colsample_bytree,
    subsample=subsample,
    reg_alpha=reg_alpha,
    reg_lambda=reg_lambda,
    min_child_weight=min_child_weight,
    gamma=gamma,
    eval_metric="rmse",
    early_stopping_rounds=50,
)

print("Entrainement avec hyperparams : ")
for param, value in xgb_model_1.get_params().items():
    if value !=None:
        print(f"{param}: {value}")
print("\n")

                                                 
# =============================================================================
# Learning model + prediction
# =============================================================================

timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M")

xgb_model_1.fit(
    X_train, y_train,
    eval_set=[(X_train, y_train), (X_test, y_test)],
    verbose=False)
evals_result_1 = xgb_model_1.evals_result()
train_rmse_1 = evals_result_1["validation_0"]["rmse"]
test_rmse_1 = evals_result_1["validation_1"]["rmse"]

plt.figure(figsize=(8,5))
plt.plot(train_rmse_1, label="Train RMSE")
plt.plot(test_rmse_1, label="Test RMSE")
plt.xlabel("n estimators")
plt.ylabel("RMSE")
plt.legend()
plt.grid(True)
plt.tight_layout()
fig_name = f"3_evolution_rmse_train_and_test_{choix_div}"
plt.savefig(f"{fig_dir}/{fig_name}", dpi=300, bbox_inches="tight")  
plt.show()

# xgb_model_1.save_model(f"/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/CODE/models_XGboost_trains/model_hz-{harmonisation_dataset}_on-{choix_div}_{timestamp}.json")

# Prediction ================================
y_pred_1 = xgb_model_1.predict(X_test)
y_train_pred_1 = xgb_model_1.predict(X_train)
        
# =============================================================================
# Evaluation
# =============================================================================

print("Evaluation metrics : ")

mse = mean_squared_error(y_test, y_pred_1)
mae = mean_absolute_error(y_test, y_pred_1)
rmse = np.sqrt(mse)
r2_test = r2_score(y_test, y_pred_1)
r2_train = r2_score(y_train, y_train_pred_1)

n_test = len(y_test)
p = X_test.shape[1]
r2_adj_test = 1 - (1 - r2_test) * (n_test - 1) / (n_test - p - 1)

n_train = len(y_train)
r2_adj_train = 1 - (1 - r2_train) * (n_train - 1) / (n_train - p - 1)
r2_adj_GAM = fit_GAM["r2_adj"][0]

print(f"MSE : {mse:.4f}") 
print(f"MAE : {mae:.4f}") 
print(f"RMSE : {rmse:.4f}")
print(f"r2 train : {r2_train:.4f}")
print(f"r2 test : {r2_test:.4f}") 
print("\n")
print(f"R2 ajusté test : {r2_adj_test:.4f}")
print(f"R2 ajusté train : {r2_adj_train:.4f}")
print(f"r2 adj GAM : {r2_adj_GAM:.4f}") 

print("\n")

# =============================================================================
# Importance des variables 
# =============================================================================

# Avec Score F = fréquence d’utilisation d’une variable comme critère de séparation dans les arbres de décision du model

fig, ax = plt.subplots(figsize=(12, 8))  
xgb.plot_importance(xgb_model_1, title="model 1", ax=ax)

fig_name = "4_Importance_VA_for_separation_tree"
fig.savefig(f"{fig_dir}/{fig_name}", dpi=300, bbox_inches="tight")
plt.show()

# fi = pd.DataFrame(data=xgb_model_1.feature_importances_,
#     index=xgb_model_1.feature_names_in_, columns=["importance"]).sort_values("importance")

# =============================================================================
# SHAP (SHapley Additive Explanations) 
# =============================================================================

# Avec valeur SHAP = indique dans quelle mesure une variable influence la prédiction du modèle. (Shap positif => influence fit positif). couleur = valeut de la variable

explainer_1 = shap.Explainer(xgb_model_1, X_train)
shap_values_1 = explainer_1(X_test)  
shap.summary_plot(shap_values_1, X_test,show=False)

# save
fig = plt.gcf()
fig_name = f"5_SHAP_{choix_div}.png"
fig.savefig(f"{fig_dir}/{fig_name}", dpi=300, bbox_inches="tight")
plt.close(fig)
shap.summary_plot(shap_values_1, X_test,show=True) # pour le réafficher

# =============================================================================
# Relation X-y
# =============================================================================

# -------- Binning sur XGBoost prédictions --------
df_XGBoost = pd.DataFrame({
    'prev_Mixo': X_test["prev_Mixo"],
    'y_pred': y_pred_1
})
df_XGBoost['bin'] = pd.cut(df_XGBoost['prev_Mixo'], bins=20)

agg_XGBoost = df_XGBoost.groupby('bin').agg(
    x_mean=('prev_Mixo', 'mean'),
    y_mean=('y_pred', 'mean'),
    y_sem=('y_pred', 'sem')
).dropna()

agg_XGBoost['ci95_hi'] = agg_XGBoost['y_mean'] + 1.96 * agg_XGBoost['y_sem']
agg_XGBoost['ci95_lo'] = agg_XGBoost['y_mean'] - 1.96 * agg_XGBoost['y_sem']

# -------- Binning sur fit GAM --------
df_gam = pd.DataFrame({
    'prev_Mixo': fit_GAM["prev_Mixo"],
    'y_fit': fit_GAM["fit"]
})
df_gam['bin'] = pd.cut(df_gam['prev_Mixo'], bins=20)

agg_gam = df_gam.groupby('bin').agg(
    x_mean=('prev_Mixo', 'mean'),
    y_mean=('y_fit', 'mean'),
    y_sem=('y_fit', 'sem')
).dropna()

agg_gam['ci95_hi'] = agg_gam['y_mean'] + 1.96 * agg_gam['y_sem']
agg_gam['ci95_lo'] = agg_gam['y_mean'] - 1.96 * agg_gam['y_sem']

# -------- Tracé --------
plt.figure(figsize=(10,8))
plt.scatter(X_train["prev_Mixo"], y_train, label="Valeurs train", color="grey", alpha=0.1, s=15)
plt.scatter(X_test["prev_Mixo"], y_test, label="Valeurs test", color="blue", alpha=0.15, s=15)
plt.scatter(X_test["prev_Mixo"], y_pred_1, color="red", alpha=0.3, s=20)
plt.scatter(fit_GAM["prev_Mixo"], fit_GAM["fit"], color="green", alpha=0.3, s=20)

# MEAN + IC XGBoost
plt.plot(agg_XGBoost['x_mean'], agg_XGBoost['y_mean'], color='red', label="Fit XGBoost", alpha=0.8)
plt.fill_between(agg_XGBoost['x_mean'], agg_XGBoost['ci95_lo'], agg_XGBoost['ci95_hi'], color='red', alpha=0.1, label="IC 95 %% XGBoost")

# MEAN + IC GAM
plt.plot(agg_gam['x_mean'], agg_gam['y_mean'], color='green', label="Fit GAM", alpha=0.8)
plt.fill_between(agg_gam['x_mean'], agg_gam['ci95_lo'], agg_gam['ci95_hi'], color='green', alpha=0.1, label="IC 95 %% GAM")

# -------- Mise en forme --------
plt.legend()
plt.xlabel("Prévalence Mixotrophs")
plt.ylabel(f"Diversité : {choix_div}")
plt.title(
    f"R2 test : {r2_score(y_test, y_pred_1):.4f} // R2 train : {r2_score(y_train, y_train_pred_1):.4f}\n"
    f"R2 adj test : {r2_adj_test:.4f} // R2 adj train : {r2_adj_train:.4f} // R2 adj test GAM : {r2_adj_GAM:.4f}",
    fontsize=10
)

fig_name = f"6_prediction_model_{choix_div}"
plt.savefig(f"{fig_dir}/{fig_name}", dpi=300, bbox_inches="tight")  
plt.show()

# =============================================================================
# Analyse résidus
# =============================================================================

residus_1 = y_test - y_pred_1

fig = plt.figure(figsize=(8, 10))

ax1 = fig.add_subplot(3, 1, 1)
ax1.scatter(y_pred_1, residus_1, alpha=0.6)
ax1.axhline(y=0, color='r', linestyle='-')
ax1.set_xlabel("y_pred")
ax1.set_ylabel("residus")
ax1.set_title("Résidus Model XGboost (voir résidus GAM sur R)")

ax2 = fig.add_subplot(3, 1, 2)
ax2.hist(residus_1, bins=50, color="grey", alpha=0.5)

ax3 = fig.add_subplot(3, 1, 3)
sm.qqplot(residus_1, line='s', ax=ax3)

plt.tight_layout()
fig_name = f"7_residus_{choix_div}"
plt.savefig(f"{fig_dir}/{fig_name}", dpi=300, bbox_inches="tight")  
plt.show()

sys.exit()























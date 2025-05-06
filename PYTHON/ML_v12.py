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
import pandas as pd
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
compare_model = params["compare_model"]
choix_div = params["choix_div"]
learning_sep_strat = params["learning_sep_strat"]
refaire_optimisation = params["refaire_optimisation"]
input_rd = params["input_rd"]
input_op = params["input_op"]
test_size = params["test_size"]
reuse_model = params["reuse_model"]
model_file = params["model_file"]
harmonisation_dataset = params["harmonisation_dataset"]

learning_rate = params["learning_rate"]
max_depth = params["max_depth"]
n_estimators = params["n_estimators"]
colsample_bytree = params["colsample_bytree"]
subsample = params["subsample"]
reg_alpha = params["reg_alpha"]
reg_lambda = params["reg_lambda"]
min_child_weight = params["min_child_weight"]
gamma = params["gamma"]

fig_dir = f"/Users/renaudsrr/Desktop/STAGE_MTL/MODELISATION/figures/{choix_div}"
os.makedirs(fig_dir, exist_ok=True)

if compare_model == "Y" and refaire_optimisation =="N":
    all_optimisation = pd.read_csv(f"/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/Data/new_csv/optimisation_hyperparms/all_optimisation_{choix_div}.csv", index_col=0)

if reuse_model == "Y" and model_file:
    xgb_model_load = XGBRegressor()
    xgb_model_load.load_model(model_file)
    print(f"Modèle chargé depuis : {model_file}")


df_instructions = pd.DataFrame([{"choix_dataset": choix_dataset, "choix_div": choix_div}])
df_instructions.to_csv("/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/Data/Instructions/df_instructions.csv", index=False)
df_instructions.to_csv("/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/Data/Instructions/df_instructions.csv")

# =============================================================================
# Run script R 
# =============================================================================

script_CARBBAS = "/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/CODE/exploration_data_CARBBAS.R"
script_LP_NLA = "/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/CODE/exploration_data_LP&NLA.R"
script_harmz = "/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/CODE/harmonisation_data.R"
script_LP_NLA = "/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/CODE/GAM_model.R"
script_IISD_ELA = "/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/CODE/exploration_data_IISD_ELA.R"

print("Run GAM model...","\n")
subprocess.run(["Rscript", script_LP_NLA], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

file_fit_GAM = f"/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/Data/new_csv/data_python/prediction_fit_{choix_div}.csv"
fit_GAM = pd.read_csv(file_fit_GAM)

# =============================================================================
# Data
# =============================================================================

if harmonisation_dataset == "N":
    if choix_dataset == "C":

        print("Run script CARBBAS...")
        subprocess.run(["Rscript", script_CARBBAS], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)   

        data_python = pd.read_csv("/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/Data/new_csv/data_python/data_python_CARBBAS.csv",index_col=0)
        VA_expl = data_python.drop(columns=["region","lake","type_prof","rich_genus","H","one_minus_D","J"])

        # Traitement des NaN (Provisioire) ============================================

        VA_expl["biomass_tot_zoo"] = VA_expl["biomass_tot_zoo"].fillna(VA_expl["biomass_tot_zoo"].median())
        VA_expl["sech.m"] = VA_expl["sech.m"].fillna(VA_expl["sech.m"].median())
        VA_expl.to_csv("/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/Data/new_csv/data_python/VA_expl_python.csv")


    if choix_dataset == "LN":

        print("Run script LP_NLA...")
        subprocess.run(["Rscript", script_LP_NLA], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

        data_python = pd.read_csv("/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/Data/new_csv/data_python/data_python_LP_NLA.csv",index_col=1)

        if learning_sep_strat in ["mix", "strat"]:

            color_mixed = "green"
            color_strat = "blue"

            # Figure
            plt.figure(figsize=(10, 8))
            plt.scatter(
                data_python.loc[data_python["Stratification"] == "mixed", "prev_Mixo"],
                data_python.loc[data_python["Stratification"] == "mixed", choix_div],
                color=color_mixed, alpha=0.7, label="mixed"
            )
            plt.scatter(
                data_python.loc[data_python["Stratification"] == "stratified", "prev_Mixo"],
                data_python.loc[data_python["Stratification"] == "stratified", choix_div],
                color=color_strat, alpha=0.7, label="stratified"
            )
            plt.xlabel("Prévalence Mixotrophs")
            plt.ylabel(f"Diversité : {choix_div}")
            plt.legend()
            plt.tight_layout()
            fig_name = f"1_sep_mixed_stratified_{choix_div}"
            plt.savefig(f"{fig_dir}/{fig_name}", dpi=300, bbox_inches="tight")
            plt.show()

            strat_dict = {"mix": "mixed", "strat": "stratified"}
            data_python = data_python[data_python["Stratification"] == strat_dict[learning_sep_strat]]

        # VA_expl final : 
        VA_expl = data_python.drop(columns=["Survey","Stratification","Ecoregion","nb_genus_mixo","rich_genus","H","one_minus_D","J","biovol_tot","biovol_mixo"])

    if choix_dataset =="ELA":
        print("Run script ELA...")
        subprocess.run(["Rscript", script_IISD_ELA], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)   
        data_python = pd.read_csv("/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/Data/new_csv/data_python/data_python_IISD_ELA.csv",index_col=0)
        VA_expl = data_python.drop(columns=["rich_genus","H","one_minus_D","J"])

else:
    print("Run script harmonisation...","\n")
    subprocess.run(["Rscript", script_harmz], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    if choix_dataset == "C":
        print("Run script CARBBAS...")
        subprocess.run(["Rscript", script_CARBBAS], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)   
        data_python = pd.read_csv("/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/Data/new_csv/data_python/data_python_CARBBAS_harmz.csv",index_col=1)
        VA_expl = data_python.drop(columns=["region","Stratification","rich_genus","H","one_minus_D","J"])

    if choix_dataset == "LN":
        print("Run script LP_NLA...")
        subprocess.run(["Rscript", script_LP_NLA], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        data_python = pd.read_csv("/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/Data/new_csv/data_python/data_phyton_LP_NLA_harmz.csv",index_col=1)
        VA_expl = data_python.drop(columns=["region","Stratification","rich_genus","H","one_minus_D","J"])

    if choix_dataset =="ELA":
        print("pas encore d'harmonisation pour ce dataset")
        sys.exit()

# =============================================================================
# Extraction data pour ML
# =============================================================================

y = data_python[choix_div] # Variable réponse
# Séparation avec jeu learning et test, ajout random_state = 42 si besoin pour cadrer résultats predi
X_train, X_test, y_train, y_test = train_test_split(VA_expl, y, test_size=test_size,random_state = 42) 

print("var train : ", np.var(y_train), " // var test : ", np.var(y_test))
print("\n")

# =============================================================================
# Exploration data
# =============================================================================

if learning_sep_strat == "mix" :
    color_fig = color_mixed
elif learning_sep_strat == "strat" :
    color_fig = color_strat
else : 
    color_fig = "black"


plt.figure(figsize=(10, 8))

plt.scatter(X_train["prev_Mixo"],y_train,color=color_fig,alpha=0.7,label="train")
plt.scatter(X_test["prev_Mixo"],y_test,color=color_fig,alpha=0.3,label="test")
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

if compare_model == "Y" and refaire_optimisation == "Y":

    # Randomise ==================================================================
    if int(input_rd) > 0:
        print("opti Rd en cours...")


        param_dist = {
            'n_estimators': [50, 100, 300, 500, 1000, 2000,2500],
            'max_depth': [3, 6, 9, 12, 15, 20,25,30],
            'learning_rate': uniform(0.001, 0.3),  
            'subsample': uniform(0.5, 0.5),  
            'colsample_bytree': uniform(0.5, 0.5),
            'reg_lambda': uniform(0, 10),  
            'reg_alpha': uniform(0, 5),  
            'gamma' : uniform(0,0.5)
        }

        xgb_model_Rd = XGBRegressor(random_state=42)
        random_search = RandomizedSearchCV(xgb_model_Rd, param_dist, cv=5, scoring='r2', n_jobs=-1, n_iter=input_rd)
        random_search.fit(X_train, y_train)

        hyperparams_Rd = pd.DataFrame.from_dict(random_search.best_params_, orient="index", columns=["Opti_Rd"])
        ordre_hyperparams = ["learning_rate", "max_depth", "n_estimators", "colsample_bytree", "subsample", "reg_alpha", "reg_lambda", "gamma"]
        all_optimisation = hyperparams_Rd.reindex(ordre_hyperparams)

    # Bayesien Optuna ==================================================================
    if int(input_op)>0:

        print("opti Op en cours...")

        def objective(trial):
            params = {
                'n_estimators': trial.suggest_int('n_estimators', 50, 2000),
                'max_depth': trial.suggest_int('max_depth', 3, 10),
                'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.3, log=True),
                'subsample': trial.suggest_float('subsample', 0.5, 1.0),
                'colsample_bytree': trial.suggest_float('colsample_bytree', 0.5, 1.0),
                'reg_lambda': trial.suggest_float('reg_lambda', 0.0, 1.0),
                'reg_alpha': trial.suggest_float('reg_alpha', 0.0, 1.0),
                'gamma' : trial.suggest_float('gamma',0.0,0.5)
            }
            
            xgb_model_Op = XGBRegressor(**params)
            
            score = cross_val_score(xgb_model_Op, X_train, y_train, cv=5, scoring="r2",n_jobs=-1)
            return score.mean()

        study = optuna.create_study(direction='maximize')
        study.optimize(objective, n_trials=input_op)

        hyperparams_Op = pd.DataFrame.from_dict(study.best_params, orient="index", columns=["Opti_Op"])
        hyperparams_Rd = pd.DataFrame(np.zeros((len(study.best_params), 1)), columns=["Opti_Rd"])

        if int(input_rd) > 0:
            hyperparams_Rd = pd.DataFrame.from_dict(random_search.best_params_, orient="index", columns=["Opti_Rd"])
            all_optimisation = pd.merge(hyperparams_Rd, hyperparams_Op, left_index=True, right_index=True, how="outer")
        else:
            all_optimisation = hyperparams_Op

        ordre_hyperparams = ["learning_rate", "max_depth", "n_estimators", "colsample_bytree", "subsample", "reg_alpha", "reg_lambda", "gamma"]
        all_optimisation = all_optimisation.reindex(ordre_hyperparams)

        all_optimisation.to_csv(f"/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/Data/new_csv/optimisation_hyperparms/all_optimisation_{choix_div}.csv")

        print("\n")
        print(f"Optimisation hyperparams pour y = {choix_div}")
        print(all_optimisation)
        print(f"Meilleur R² Optuna : {study.best_value:.4f}")
        print("\n")

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

if compare_model=="N":
    print("Entrainement avec hyperparams : ")
    for param, value in xgb_model_1.get_params().items():
        if value !=None:
            print(f"{param}: {value}")
    print("\n")

if compare_model =="Y": 

    if int(input_rd)>0:

        xgb_model_2 = XGBRegressor(
            random_state = 42,        
            # grow_policy = "depthwise",
            # importance_type = "gain",

            learning_rate = all_optimisation.loc["learning_rate", "Opti_Rd"],    
            max_depth = int(all_optimisation.loc["max_depth", "Opti_Rd"]),          
            n_estimators = int(all_optimisation.loc["n_estimators", "Opti_Rd"]),    

            colsample_bytree = all_optimisation.loc["colsample_bytree", "Opti_Rd"],  
            subsample = all_optimisation.loc["subsample", "Opti_Rd"],         

            reg_alpha = all_optimisation.loc["reg_alpha", "Opti_Rd"],  
            reg_lambda = all_optimisation.loc["reg_lambda", "Opti_Rd"],
            
            gamma = all_optimisation.loc["gamma", "Opti_Rd"],  

            early_stopping_rounds=50,
            )
    
    if int(input_op)>0:

        xgb_model_3 = XGBRegressor(
            random_state = 42,        
            # grow_policy = "depthwise",
            # importance_type = "gain",
            
            learning_rate = all_optimisation.loc["learning_rate", "Opti_Op"],    
            max_depth =  int(all_optimisation.loc["max_depth", "Opti_Op"]),          
            n_estimators =  int(all_optimisation.loc["n_estimators", "Opti_Op"]),    

            colsample_bytree =  all_optimisation.loc["colsample_bytree", "Opti_Op"],  
            subsample =  all_optimisation.loc["subsample", "Opti_Op"],         

            reg_alpha =  all_optimisation.loc["reg_alpha", "Opti_Op"],  
            reg_lambda =  all_optimisation.loc["reg_lambda", "Opti_Op"],
            
            gamma = all_optimisation.loc["gamma", "Opti_Op"],  

            early_stopping_rounds=50,
            )
                                                 
# =============================================================================
# Learning model + prediction
# =============================================================================

timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M")

if compare_model == "Y":
    xgb_model_1.fit(
        X_train, y_train,
        eval_set=[(X_train, y_train), (X_test, y_test)],
        verbose=False)
    evals_result_1 = xgb_model_1.evals_result()
    train_rmse_1 = evals_result_1["validation_0"]["rmse"]
    test_rmse_1 = evals_result_1["validation_1"]["rmse"]
    all_rmse = train_rmse_1 + test_rmse_1

    # xgb_model_1.save_model(f"models_trains/hz_{harmonisation_dataset}_{choix_div}_{timestamp}.json")

    if int(input_rd)>0:
        xgb_model_2.fit(    
            X_train, y_train,
            eval_set=[(X_train, y_train), (X_test, y_test)],
            verbose=False)
        evals_result_2 = xgb_model_2.evals_result()
        train_rmse_2 = evals_result_2["validation_0"]["rmse"]
        test_rmse_2 = evals_result_2["validation_1"]["rmse"]
        all_rmse += train_rmse_2 + test_rmse_2

        # xgb_model_2.save_model(f"models_trains/xgb_model_2_{choix_div}_{timestamp}.json")
    
    if int(input_op)>0:
        xgb_model_3.fit(
            X_train, y_train,
            eval_set=[(X_train, y_train), (X_test, y_test)],
            verbose=False)
        evals_result_3 = xgb_model_3.evals_result()
        train_rmse_3 = evals_result_3["validation_0"]["rmse"]
        test_rmse_3 = evals_result_3["validation_1"]["rmse"]
        all_rmse += train_rmse_3 + test_rmse_3
        
        # xgb_model_3.save_model(f"models_trains/xgb_model_3_{choix_div}_{timestamp}.json")

    min_rmse = min(all_rmse)
    max_rmse = max(all_rmse)

    plt.figure(figsize=(11,8))

    plt.subplot(1,1+int(input_rd)+int(input_op),1)
    plt.plot(train_rmse_1, label="Train RMSE")
    plt.plot(test_rmse_1, label="Test RMSE")
    plt.xlabel("n estimators")
    plt.ylabel("RMSE")
    plt.title("Model 1")
    plt.ylim(min_rmse, max_rmse)
    plt.legend()
    plt.grid(True)

    if int(input_rd) > 0:
        plt.subplot(1,1+int(input_rd)+int(input_op),2)
        plt.plot(train_rmse_2, label="Train RMSE")
        plt.plot(test_rmse_2, label="Test RMSE")
        plt.xlabel("n estimators")
        plt.ylabel("RMSE")
        plt.title("Model 2")
        plt.ylim(min_rmse, max_rmse)
        plt.legend()
        plt.grid(True)

        # Prediction ================================
        y_pred_2 = xgb_model_2.predict(X_test)
        y_train_pred_2 = xgb_model_2.predict(X_train)

    if int(input_op) > 0:
        plt.subplot(1,1+int(input_rd)+int(input_op),2 if int(input_rd) == 0 else 3)
        plt.plot(train_rmse_3, label="Train RMSE")
        plt.plot(test_rmse_3, label="Test RMSE")
        plt.xlabel("n estimators")
        plt.ylabel("RMSE")
        plt.title("Model 3")
        plt.ylim(min_rmse, max_rmse)
        plt.legend()
        plt.grid(True)

        # Prediction ================================
        y_pred_3 = xgb_model_3.predict(X_test)
        y_train_pred_3 = xgb_model_3.predict(X_train)

    # Prediction ================================
    y_pred_1 = xgb_model_1.predict(X_test)
    y_train_pred_1 = xgb_model_1.predict(X_train)

    plt.tight_layout()
    plt.show()

else:

    if reuse_model =="N":

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

        xgb_model_1.save_model(f"/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/models_trains/model_hz-{harmonisation_dataset}_on-{choix_div}_{timestamp}.json")

        # Prediction ================================
        y_pred_1 = xgb_model_1.predict(X_test)
        y_train_pred_1 = xgb_model_1.predict(X_train)

    else:

        # Prediction ================================
        y_pred_1 = xgb_model_load.predict(X_test)
        y_train_pred_1 = xgb_model_load.predict(X_train)

# =============================================================================
# Evaluation
# =============================================================================

print("Evaluation metrics : ")

if compare_model =="Y":

    if int(input_rd) == 0:
        y_pred_2 = np.zeros(y_pred_1.shape[0])
        y_train_pred_2 = np.zeros(y_train_pred_1.shape[0])
    if int(input_op) == 0: 
        y_pred_3 = np.zeros(y_pred_1.shape[0])
        y_train_pred_3 = np.zeros(y_train_pred_1.shape[0])

    metrics_eval = {
        "Model": ["Model_1", "Model_2", "Model_3"],
        "MSE": [
            mean_squared_error(y_test, y_pred_1),
            mean_squared_error(y_test, y_pred_2),
            mean_squared_error(y_test, y_pred_3)
        ],
        "MAE": [
            mean_absolute_error(y_test, y_pred_1),
            mean_absolute_error(y_test, y_pred_2),
            mean_absolute_error(y_test, y_pred_3)
        ],
        "RMSE": [
            np.sqrt(mean_squared_error(y_test, y_pred_1)),
            np.sqrt(mean_squared_error(y_test, y_pred_2)),
            np.sqrt(mean_squared_error(y_test, y_pred_3))
        ],
        "R2 test": [
            r2_score(y_test, y_pred_1),
            r2_score(y_test, y_pred_2),
            r2_score(y_test, y_pred_3)
        ],
        "R2 train": [
            r2_score(y_train, y_train_pred_1),
            r2_score(y_train, y_train_pred_2),
            r2_score(y_train, y_train_pred_3)
        ]
    }
    metric_eval_df = pd.DataFrame(metrics_eval)
    print(metric_eval_df)   

else : 
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
if compare_model =="Y":
    xgb.plot_importance(xgb_model_1,title = "model 1")

    if int(input_rd)>0:
        xgb.plot_importance(xgb_model_2,title = "model 2")
    if int(input_op)>0:
        xgb.plot_importance(xgb_model_3,title = "model 3")

    plt.show()

else :
    if reuse_model == "N":

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

if compare_model =="Y":
    explainer_1 = shap.Explainer(xgb_model_1, X_train)
    shap_values_1 = explainer_1(X_test)
    shap.summary_plot(shap_values_1, X_test)

    if int(input_rd)>0:    
        explainer_2 = shap.Explainer(xgb_model_2, X_train)
        shap_values_2 = explainer_2(X_test)
        shap.summary_plot(shap_values_2, X_test)

    if int(input_op)>0:    
        explainer_3 = shap.Explainer(xgb_model_3, X_train)
        shap_values_3 = explainer_3(X_test)
        shap.summary_plot(shap_values_3, X_test)
else : 
    if reuse_model == "N":

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

sorted_indices = np.argsort(X_test["prev_Mixo"])
X_sorted = X_test["prev_Mixo"].iloc[sorted_indices]

y_pred_sorted_1 = y_pred_1[sorted_indices]

if compare_model =="Y":
    y_pred_sorted_2 = y_pred_2[sorted_indices]
    y_pred_sorted_3 = y_pred_3[sorted_indices]

plt.figure(figsize=(10,8))

if compare_model =="Y":
    plt.subplot(1,3,1)

plt.scatter(X_train["prev_Mixo"], y_train, label="Valeurs train", color="grey", alpha=0.2)
plt.scatter(X_test["prev_Mixo"], y_test, label="Valeurs test réelles", color="blue", alpha=0.4)
plt.scatter(X_test["prev_Mixo"], y_pred_1, label="Valeurs prédites model 1", color="red", alpha=0.45)
plt.step(X_sorted, y_pred_sorted_1, where="post", color="red", linestyle="-", linewidth=1, alpha=0.5)

if choix_dataset!="ELA":
    plt.scatter(fit_GAM["prev_Mixo"],fit_GAM["fit"], label="fit GAM", color="green", alpha=0.4)

plt.legend()
plt.xlabel("Prévalence Mixotrophs")
plt.ylabel(f"Diversite : {choix_div}")

if choix_dataset !="ELA":
    plt.title(
        f"R2 test : {r2_score(y_test, y_pred_1):.4f} // R2 train : {r2_score(y_train, y_train_pred_1):.4f}\n"
        f"R2 adj test : {r2_adj_test:.4f} // R2 adj train : {r2_adj_train:.4f} // R2 adj test GAM : {r2_adj_GAM:.4f}",
        fontsize=10
    )

if compare_model =="Y":
    plt.subplot(1,3,2)
    plt.scatter(X_train["prev_Mixo"], y_train, label="Valeurs train", color="grey", alpha=0.2)
    plt.scatter(X_test["prev_Mixo"], y_test, label="Valeurs test réelles", color="blue", alpha=0.4)
    plt.scatter(X_test["prev_Mixo"], y_pred_2, label="Valeurs prédites model 2", color="orange", alpha=0.45)
    plt.step(X_sorted, y_pred_sorted_2, where="post", color="orange", linestyle="-", linewidth=1, alpha=0.5)
    plt.legend()
    plt.xlabel("Prévalence Mixotrophs")
    plt.title(f"R2 test : {r2_score(y_test, y_pred_2):.4f} // R2 train : {r2_score(y_train, y_train_pred_2):.4f}", fontsize=10)

    plt.subplot(1,3,3)
    plt.scatter(X_train["prev_Mixo"], y_train, label="Valeurs train", color="grey", alpha=0.2)
    plt.scatter(X_test["prev_Mixo"], y_test, label="Valeurs test réelles", color="blue", alpha=0.4)
    plt.scatter(X_test["prev_Mixo"], y_pred_3, label="Valeurs prédites model 3", color="purple", alpha=0.45)
    plt.step(X_sorted, y_pred_sorted_3, where="post", color="purple", linestyle="-", linewidth=1, alpha=0.5)
    plt.legend()
    plt.xlabel("Prévalence Mixotrophs")
    plt.title(f"R2 test : {r2_score(y_test, y_pred_3):.4f} // R2 train : {r2_score(y_train, y_train_pred_3):.4f}", fontsize=10)

    plt.tight_layout()

fig_name = f"6_prediction_model_{choix_div}"
plt.savefig(f"{fig_dir}/{fig_name}", dpi=300, bbox_inches="tight")  
plt.show()


# =============================================================================
# Analyse résidus
# =============================================================================

if compare_model =="Y":
    residus_1 = y_test - y_pred_1
    residus_2 = y_test - y_pred_2
    residus_3 = y_test - y_pred_3

    plt.subplot(3,2,1)
    plt.scatter(y_pred_1, residus_1, alpha=0.6)
    plt.axhline(y=0, color='r', linestyle='-', alpha=0.2)
    plt.xlabel("y_pred")
    plt.ylabel("residus_1")
    plt.subplot(3,2,2)
    plt.hist(residus_1, bins=30,color="grey",alpha=0.5)

    plt.subplot(3,2,3)
    plt.scatter(y_pred_2, residus_2, alpha=0.6)
    plt.axhline(y=0, color='r', linestyle='-', alpha=0.2)
    plt.xlabel("y_pred")
    plt.ylabel("residus_2")
    plt.subplot(3,2,4)
    plt.hist(residus_2, bins=30,color="grey",alpha=0.5)

    plt.subplot(3,2,5)
    plt.scatter(y_pred_3, residus_3, alpha=0.6)
    plt.axhline(y=0, color='r', linestyle='-', alpha=0.2)
    plt.xlabel("y_pred")    
    plt.ylabel("residus_3")
    plt.subplot(3,2,6)
    plt.hist(residus_3, bins=30, color="grey",alpha=0.5)
    
    plt.tight_layout()
    plt.show()

    print("NA dans : ")
    print(X_train.isna().sum())         # combien de NaN par variable


else:
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























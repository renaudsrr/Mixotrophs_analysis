"""
Script building a model [...]

Author : SERRE Renaud
Creation 12.03.2025

History of modification
12.03.2025 : Creation of the script. Import data 
13.03.2025 : Some test of RF
18.03.2025 : Contnuation and basics stuff of XGboost, SHAP

"""

# =============================================================================
# Packages
# =============================================================================

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from datetime import datetime

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

# =======================================================================================================================
# =======================================================================================================================
# PARAM SCRIPT 
# =======================================================================================================================
# =======================================================================================================================

# fonction pour reuse_model
def choice_model():
    dossier = "/Users/renaudsrr/Desktop/STAGE_MTL/Scripts/models_trains"
    fichiers = [f for f in os.listdir(dossier) if os.path.isfile(os.path.join(dossier, f))]

    print("Fichiers disponibles :")
    for i, fichier in enumerate(fichiers):
        print(f"{i + 1}. {fichier}")
    
    while True:
        choix = int(input("Num fichier ? : ")) - 1
        if 0 <= choix < len(fichiers):
            fichier_choisi = fichiers[choix]
            return os.path.join(dossier, fichier_choisi)



# Choix de paramétrage de modélisation
choix_dataset = input("(Vide = LN) // Choix dataset utilisé : (LN/C)") or "LN"

compare_model = input("(Vide = N) // Comparaison model simple + 2 model optimisés ? (Y/N)") or "N"

choix_div = input("(Vide = H) // Choix indice de diversite ? (rich_genus, H, 1-D, J) ") or "H"

if choix_dataset =="LN":
    learning_sep_strat = input("(Vide = N) // Séparer le jeu d'entrainement en fct stratification ? (N/mix/strat)") or "N"
else:
    learning_sep_strat = "N"

# reuse_model = input("(Vide = N) // Réutiliser un model performant déjà entrainé ? (Y/N)") or "N"
# if reuse_model =="Y":
#     model_path = choice_model()
#     xgb_model_load = xgb.XGBRegressor()
#     xgb_model_load.load_model(model_path)
#     print(f"Modèle chargé depuis : {model_path}")
reuse_model = "N"

if compare_model =="Y":
    refaire_optimisation = input("(Vide = N) // Refaire l'optimisation hyperparams ? (Y/N) ") or "N"
    if refaire_optimisation =="Y":
        input_rd = input("(Vide = 1) // Nb essais optimisation hyperparam Randomised : ")
        n_trials_opt_Rd = int(input_rd) if input_rd.strip() != "" else 1

        input_op = input("(Vide = 1) // Nb essais optimisation hyperparam optuna : ")
        n_trials_opt_Op = int(input_op) if input_op.strip() != "" else 1
    else:
        all_optimisation = pd.read_csv(f"Data/new_csv/optimisation_hyperparms/all_optimisation_{choix_div}.csv", index_col=0)
        input_rd = int("Opti_Rd" in all_optimisation.columns)
        input_op = int("Opti_Op" in all_optimisation.columns)

if compare_model == "Y" and refaire_optimisation =="N":
    all_optimisation = pd.read_csv(f"Data/new_csv/optimisation_hyperparms/all_optimisation_{choix_div}.csv", index_col=0)

# % données utilisées pour X_test et y_test
test_size = 0.2

# pour enregistrement img figures
fig_dir = f"figures/{choix_div}"
os.makedirs(fig_dir, exist_ok=True)

print("\n")

# =============================================================================
# Data
# =============================================================================

if choix_dataset == "C":
    data_python = pd.read_csv("Data/new_csv/data_python/data_python_CARBBAS.csv",index_col=0)
    VA_expl = data_python.drop(columns=["region","lake","type_prof","rich_genus","H","1-D","J"])

    # Traitement des NaN (Provisioire) ============================================

    VA_expl["biomass_tot_zoo"] = VA_expl["biomass_tot_zoo"].fillna(VA_expl["biomass_tot_zoo"].median())
    VA_expl["sech.m"] = VA_expl["sech.m"].fillna(VA_expl["sech.m"].median())
    VA_expl.to_csv("Data/new_csv/data_python/VA_expl_python.csv")


if choix_dataset == "LN":
    data_python = pd.read_csv("Data/new_csv/data_python/data_python_LP_NLA.csv",index_col=1)

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
        fig_name = "1_sep_mixed_stratified"
        plt.savefig(f"{fig_dir}/{fig_name}", dpi=300, bbox_inches="tight")
        plt.show()

        strat_dict = {"mix": "mixed", "strat": "stratified"}
        data_python = data_python[data_python["Stratification"] == strat_dict[learning_sep_strat]]

    # VA_expl final : 
    VA_expl = data_python.drop(columns=["Survey","Stratification","Ecoregion","nb_genus_mixo","rich_genus","H","1-D","J"])

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
        random_search = RandomizedSearchCV(xgb_model_Rd, param_dist, cv=5, scoring='r2', n_jobs=-1, n_iter=n_trials_opt_Rd)
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
        study.optimize(objective, n_trials=n_trials_opt_Op)

        hyperparams_Op = pd.DataFrame.from_dict(study.best_params, orient="index", columns=["Opti_Op"])
        hyperparams_Rd = pd.DataFrame(np.zeros((len(study.best_params), 1)), columns=["Opti_Rd"])

        if int(input_rd) > 0:
            hyperparams_Rd = pd.DataFrame.from_dict(random_search.best_params_, orient="index", columns=["Opti_Rd"])
            all_optimisation = pd.merge(hyperparams_Rd, hyperparams_Op, left_index=True, right_index=True, how="outer")
        else:
            all_optimisation = hyperparams_Op

        ordre_hyperparams = ["learning_rate", "max_depth", "n_estimators", "colsample_bytree", "subsample", "reg_alpha", "reg_lambda", "gamma"]
        all_optimisation = all_optimisation.reindex(ordre_hyperparams)

        all_optimisation.to_csv(f"Data/new_csv/optimisation_hyperparms/all_optimisation_{choix_div}.csv")

        print("\n")
        print(f"Optimisation hyperparams pour y = {choix_div}")
        print(all_optimisation)
        print(f"Meilleur R² Optuna : {study.best_value:.4f}")
        print("\n")

# =============================================================================
# Creation model
# =============================================================================

xgb_model_1 = XGBRegressor(
    random_state = 42,       # Graine aléatoire pour assurer la reproductibilité des résultats

    # grow_policy = "depthwise",   # Politique de croissance des arbres (par défaut : "depthwise", sinon "lossguide")

    # importance_type = "gain",   # Mesure l'importance des variables en fonction du gain apporté par chaque split

    learning_rate = 0.008,    # Taux d'apprentissage (plus bas = apprentissage lent mais meilleure généralisation, "shrinkage")
    max_depth = 1,           # Profondeur max des arbres (plus élevé = plus de complexité et risque d'overfitting)
    n_estimators = 2000,      # Nombre total d'arbres dans l'ensemble (plus élevé = modèle plus complexe et plus lent)

    colsample_bytree = 0.8,   # Proportion des colonnes (features) utilisées pour chaque arbre (ex: 80% des variables)
    subsample = 0.8,          # Proportion des échantillons (lignes) utilisées pour chaque arbre (ex: 80% des lignes) 
    #                          # Permet d'ajouter de la variance et réduit le risque d'overfitting

    reg_alpha = 0.5,         # Régularisation L1 (LASSO) : Encourage la réduction de certains coefficients à zéro (sparse model)
    reg_lambda = 3.0,         # Régularisation L2 (Ridge) : Réduit la magnitude des coefficients pour éviter les valeurs extrêmes

    min_child_weight = 10,    # Nombre minimal d'échantillons requis dans un nœud pour qu'il puisse être divisé 
                             # (plus haut = moins de splits, prévient l'overfitting)

    gamma = 0.3,            # Seuil minimum de réduction de la perte pour effectuer un split 
                             # (plus élevé = moins de splits inutiles)

    # objective = "reg:squarederror",   # Fonction de coût pour un problème de régression classique (MSE)

    # booster = "gbtree",      # Type de booster utilisé (gbtree = gradient boosting tree)

    eval_metric = "rmse",    # Indicateur d'évaluation utilisé (par défaut)
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

    xgb_model_1.save_model(f"models_trains/xgb_model_1_{choix_div}_{timestamp}.json")

    if int(input_rd)>0:
        xgb_model_2.fit(    
            X_train, y_train,
            eval_set=[(X_train, y_train), (X_test, y_test)],
            verbose=False)
        evals_result_2 = xgb_model_2.evals_result()
        train_rmse_2 = evals_result_2["validation_0"]["rmse"]
        test_rmse_2 = evals_result_2["validation_1"]["rmse"]
        all_rmse += train_rmse_2 + test_rmse_2

        xgb_model_2.save_model(f"models_trains/xgb_model_2_{choix_div}_{timestamp}.json")
    
    if int(input_op)>0:
        xgb_model_3.fit(
            X_train, y_train,
            eval_set=[(X_train, y_train), (X_test, y_test)],
            verbose=False)
        evals_result_3 = xgb_model_3.evals_result()
        train_rmse_3 = evals_result_3["validation_0"]["rmse"]
        test_rmse_3 = evals_result_3["validation_1"]["rmse"]
        all_rmse += train_rmse_3 + test_rmse_3
        
        xgb_model_3.save_model(f"models_trains/xgb_model_3_{choix_div}_{timestamp}.json")

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
        fig_name = "3_evolution_rmse_train_and_test"
        plt.savefig(f"{fig_dir}/{fig_name}", dpi=300, bbox_inches="tight")  
        plt.show()

        xgb_model_1.save_model(f"models_trains/xgb_model_{choix_div}_{timestamp}.json")

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

    print(f"MSE : {mse:.4f}") 
    print(f"MAE : {mae:.4f}") 
    print(f"RMSE : {rmse:.4f}")
    print(f"R2 test : {r2_test:.4f}") # part de variance expliqué par le model
    print(f"R2 train : {r2_train:.4f}")

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
    explainer_1 = shap.Explainer(xgb_model_1, X_train)
    shap_values_1 = explainer_1(X_test)  
    shap.summary_plot(shap_values_1, X_test,show=False)

    # save
    fig = plt.gcf()
    fig_name = "5_SHAP.png"
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
plt.legend()
plt.xlabel("Prévalence Mixotrophs")
plt.ylabel(f"Diversite : {choix_div}")
plt.title(f"R2 test : {r2_score(y_test, y_pred_1):.4f} // R2 train : {r2_score(y_train, y_train_pred_1):.4f}", fontsize=10)

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

fig_name = "6_prediction_model"
plt.savefig(f"{fig_dir}/{fig_name}", dpi=300, bbox_inches="tight")  
plt.show()


# =============================================================================
# Analyse résidus
# =============================================================================

if compare_model =="Y":
    residus_1 = y_test - y_pred_1
    residus_2 = y_test - y_pred_2
    residus_3 = y_test - y_pred_3

    plt.subplot(4,2,1)
    plt.scatter(y_pred_1, residus_1, alpha=0.6)
    plt.axhline(y=0, color='r', linestyle='-', alpha=0.2)
    plt.xlabel("y_pred")
    plt.ylabel("residus_1")
    plt.subplot(4,2,2)
    plt.hist(residus_1, bins=30,color="grey",alpha=0.5)

    plt.subplot(4,2,3)
    plt.scatter(y_pred_2, residus_2, alpha=0.6)
    plt.axhline(y=0, color='r', linestyle='-', alpha=0.2)
    plt.xlabel("y_pred")
    plt.ylabel("residus_2")
    plt.subplot(4,2,4)
    plt.hist(residus_2, bins=30,color="grey",alpha=0.5)

    plt.subplot(4,2,5)
    plt.scatter(y_pred_3, residus_3, alpha=0.6)
    plt.axhline(y=0, color='r', linestyle='-', alpha=0.2)
    plt.xlabel("y_pred")    
    plt.ylabel("residus_3")
    plt.subplot(4,2,6)
    plt.hist(residus_3, bins=30, color="grey",alpha=0.5)
    
    plt.subplot(4,2,7)


    X_train_const = sm.add_constant(X_train) 
    df_ols = pd.concat([pd.DataFrame(X_train_const), pd.Series(y_train, name='y')], axis=1)
    df_ols_clean = df_ols.replace([np.inf, -np.inf], np.nan).dropna()
    X_train_clean = df_ols_clean.drop(columns='y')
    y_train_clean = df_ols_clean['y']
    model = sm.OLS(y_train_clean, X_train_clean).fit()

    print("NA dans : ")
    print(X_train_const.isna().sum())         # combien de NaN par variable

    influence = OLSInfluence(model)
    cooks_d = influence.cooks_distance[0]

    plt.stem(cooks_d, markerfmt=",")
    plt.ylabel("Cook’s Distance")

    plt.tight_layout()
    plt.show()

else:
    residus_1 = y_test - y_pred_1

    plt.subplot(2,1,1)
    plt.scatter(y_pred_1, residus_1, alpha=0.6)
    plt.axhline(y=0, color='r', linestyle='-')
    plt.xlabel("y_pred")
    plt.ylabel("residus_1")

    plt.subplot(2,1,2)
    plt.hist(residus_1, bins=50,color="grey",alpha=0.5)

    plt.tight_layout()
    fig_name = "7_residus"
    plt.savefig(f"{fig_dir}/{fig_name}", dpi=300, bbox_inches="tight")  
    plt.show()





























"""
Script building a model [...]

Author : SERRE Renaud
Creation 12.03.2025

History of modification
12.03.2025 : Creation of the script. Import data 
13.03.2025 : Some test of RF
18.03.2025 : Contnuation and basics stuff of XGboost, SHAP
20.03.2025 : Huge amelioration implementation of different choice of modelisation
"""

# =============================================================================
# Packages
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

from scipy.stats import uniform,stats
from statsmodels.stats.outliers_influence import OLSInfluence

import xgboost as xgb
import shap
import optuna

from sklearn.model_selection import train_test_split, GridSearchCV, RandomizedSearchCV
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from sklearn.inspection import PartialDependenceDisplay
from sklearn.preprocessing import StandardScaler

# =======================================================================================================================
# =======================================================================================================================
# PARAM SCRIPT 
# =======================================================================================================================
# =======================================================================================================================

compare_model = input("(Vide = N) // Comparaison 3 model ? (Y/N)") or "N"

if compare_model =="Y":
    refaire_optimisation = input("(Vide = N) // Refaire l'optimisation hyperparams ? (Y/N) ") or "N"
    if refaire_optimisation =="Y":
        n_trials_opt_Rd = int(input("(Vide = 100) // Nb essais optimisation hyperparam Randomised : ")) or 100
        n_trials_opt_Op = int(input("(Vide = 80) // Nb essais optimisation hyperparam optuna : ")) or 80

choix_div = input("(Vide = H) // Choix indice de diversite ? (rich_genus, H, 1-D, J) ") or "H"

normalisation_VA_expl = input("(Vide = N) // Normaliser les VA explicatives ? (Y/N) ") or "N"

test_size = 0.3

print("\n")

# =============================================================================
# Data
# =============================================================================

data_phyto = pd.read_csv("Data/new_csv/data_phyto.csv") 
data_zoo = pd.read_csv("Data/new_csv/data_zoo.csv")

# DATA AVEC : 1ère VA Reponse = 4 métriques de diversité phyto
summary_data_phyto = pd.read_csv("Data/new_csv/summary_data_phyto.csv") 
#  DATA AVEC : 2ème VA Reponse = 4 métriques de diversité zoo
summary_data_zoo = pd.read_csv("Data/new_csv/summary_data_zoo.csv") 

# DATA AVEC : all VA Explicatives
VA_expl = pd.read_csv("Data/new_csv/VA_expl.csv")
VA_expl = VA_expl.drop(columns=["Sample","region","lake","type_prof","doy","yr"])

# =============================================================================
# Exploration data
# =============================================================================

plt.figure(figsize=(10, 8))

plt.subplot(2, 2, 1)
plt.scatter(summary_data_phyto["prev_Mixo"], summary_data_phyto["rich_genus"], alpha=0.7, color="blue")
plt.xlabel("prev_Mixo")
plt.ylabel("rich_genus")

plt.subplot(2, 2, 2)
plt.scatter(summary_data_phyto["prev_Mixo"], summary_data_phyto["H"], alpha=0.7, color="green")
plt.xlabel("prev_Mixo")
plt.ylabel("H")

plt.subplot(2, 2, 3)
plt.scatter(summary_data_phyto["prev_Mixo"], summary_data_phyto["1-D"], alpha=0.7, color="red")
plt.xlabel("prev_Mixo")
plt.ylabel("1-D")

plt.subplot(2, 2, 4)
plt.scatter(summary_data_phyto["prev_Mixo"], summary_data_phyto["J"], alpha=0.7, color="purple")
plt.xlabel("prev_Mixo")
plt.ylabel("J")

plt.tight_layout()
plt.show()

# =============================================================================
# Traitement des NaN (Provisioire)
# =============================================================================

VA_expl["biomass_tot_zoo"] = VA_expl["biomass_tot_zoo"].fillna(VA_expl["biomass_tot_zoo"].median())
VA_expl["sech.m"] = VA_expl["sech.m"].fillna(VA_expl["sech.m"].median())

VA_expl.to_csv("Data/new_csv/VA_expl_python.csv", index=False)

# =============================================================================
# Extraction data pour ML
# =============================================================================

y = summary_data_phyto[choix_div] # Variable réponse
# Séparation avec jeu learning et test, ajout random_state = 42 si besoin pour cadrer résultats predi
X_train, X_test, y_train, y_test = train_test_split(VA_expl, y, test_size=test_size) 

# =============================================================================
# Normalisation
# =============================================================================

if normalisation_VA_expl =="Y":
    scaler = StandardScaler()

    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    X_train = pd.DataFrame(X_train_scaled, columns=X_train.columns, index=X_train.index)
    X_test = pd.DataFrame(X_test_scaled, columns=X_test.columns, index=X_test.index)

# =============================================================================
# Optimisation hyperparamètres
# =============================================================================

if compare_model == "Y" and refaire_optimisation == "Y":

    # Randomise ==================================================================

    param_dist = {
        'n_estimators': [50, 100, 300, 500, 1000, 2000],
        'max_depth': [3, 6, 9, 12, 15, 20],
        'learning_rate': uniform(0.001, 0.3),  
        'subsample': uniform(0.5, 0.5),  
        'colsample_bytree': uniform(0.5, 0.5),
        'reg_lambda': uniform(0, 10),  
        'reg_alpha': uniform(0, 5),  
        'min_child_weight': [1, 3, 5, 7, 10],
        'gamma': uniform(0, 5)
    }

    xgb_model_Rd = xgb.XGBRegressor(random_state=42)
    random_search = RandomizedSearchCV(xgb_model_Rd, param_dist, cv=5, scoring='r2', n_jobs=-1, n_iter=n_trials_opt_Rd, verbose=3)
    random_search.fit(X_train, y_train)

    # Bayesien Optuna ==================================================================

    def objective(trial):
        params = {
            'n_estimators': trial.suggest_int('n_estimators', 50, 3000),
            'max_depth': trial.suggest_int('max_depth', 3, 20),
            'learning_rate': trial.suggest_float('learning_rate', 0.001, 0.3, log=True),
            'subsample': trial.suggest_float('subsample', 0.5, 1.0),
            'colsample_bytree': trial.suggest_float('colsample_bytree', 0.5, 1.0),
            'random_state': 42,
            'reg_lambda': trial.suggest_float('reg_lambda', 0.0, 10.0),
            'reg_alpha': trial.suggest_float('reg_alpha', 0.0, 5.0),
            'booster': 'gbtree',
            'objective': 'reg:squarederror'
        }
        
        xgb_model_Op = xgb.XGBRegressor(**params)
        xgb_model_Op.fit(X_train, y_train)
        preds = xgb_model_Op.predict(X_test)
        return r2_score(y_test, preds)

    study = optuna.create_study(direction='maximize')
    study.optimize(objective, n_trials=n_trials_opt_Op)

    hyperparams_Rd = pd.DataFrame.from_dict(random_search.best_params_, orient="index", columns=["Opti_Rd"])
    hyperparams_Op = pd.DataFrame.from_dict(study.best_params, orient="index", columns=["Opti_Op"])

    all_optimisation = pd.merge(hyperparams_Rd, hyperparams_Op, left_index=True, right_index=True, how="outer")
    ordre_hyperparams = ["learning_rate", "max_depth", "n_estimators", "colsample_bytree", "subsample", "reg_alpha", "reg_lambda", "min_child_weight", "gamma", "max_delta_step", "booster"]
    all_optimisation = all_optimisation.reindex(ordre_hyperparams)

    all_optimisation.to_csv(f"Data/new_csv/all_optimisation_{choix_div}.csv")

    print("\n\n")
    print(f"Optimisation hyperparams pour y = {choix_div}")
    print(all_optimisation)
    print("\n\n")

# =============================================================================
# Creation model
# =============================================================================

xgb_model_1 = xgb.XGBRegressor(
    random_state = 42,        # graine aleatoire pour reproductibilite
    # grow_policy = "depthwise",
    # importance_type = "gain",

    eta = 0.01,

    learning_rate = 0.1,    # Plus il est bas, plus le modèle apprend lentement mais en généralisant mieux => Shrinkage
    max_depth = 6,           # Profondeur max des arbres => traux haut = overfitting
    n_estimators = 800,     # nb arbres

    colsample_bytree = 0.8,  # Prend 80 % des colonnes = variables dans chaque arbre
    subsample = 0.8,         # Prend 80 % des lignes du jeu learning pour chaque arbrew => evite overfitting

    reg_alpha = 0.05,  # L1 Regularization : LASSO : Favorise des coefficients proches de zero et peut eliminer certaines variables inutiles
    reg_lambda = 0.5,  # L2 Regularization : RIDGE Reduit l'importance des coefficients pour eviter les valeurs extremes
    
    min_child_weight = 4,  # Permet de faire plus de splits, réduit le sous-ajustement
    gamma = 0.2,  # Prune légèrement les splits inutiles, mais garde une bonne flexibilité

    objective = "reg:squarederror",  # Objectif de régression classique
    booster = "gbtree",  # Utilisation de l’approche en arbres
    importance_type = "gain",  # Importance des features basée sur le gain

    eval_metric = "auc"
    )

if compare_model=="N":
    print("Entrainement avec hyperparams : ")
    for param, value in xgb_model_1.get_params().items():
        if value !=None:
            print(f"{param}: {value}")
    print("\n")

if compare_model =="Y": 

    if refaire_optimisation =="N": 
        all_optimisation = pd.read_csv(f"Data/new_csv/all_optimisation_{choix_div}.csv",index_col=0)

    xgb_model_2 = xgb.XGBRegressor(
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
        )

    xgb_model_3 = xgb.XGBRegressor(
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
        )
                                             
# =============================================================================
# Learning model
# =============================================================================

if compare_model=="Y":
    xgb_model_1.fit(X_train,y_train)
    xgb_model_2.fit(X_train,y_train)
    xgb_model_3.fit(X_train,y_train)
else:
    xgb_model_1.fit(X_train,y_train)

# =============================================================================
# Prediction
# =============================================================================

if compare_model =="Y":
    y_pred_1 = xgb_model_1.predict(X_test)
    y_pred_2 = xgb_model_2.predict(X_test)
    y_pred_3 = xgb_model_3.predict(X_test)
else : 
    y_pred_1 = xgb_model_1.predict(X_test)


# =============================================================================
# Evaluation
# =============================================================================

print("Evaluation metrics : ")

if compare_model =="Y":
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
        "R2": [
            r2_score(y_test, y_pred_1),
            r2_score(y_test, y_pred_2),
            r2_score(y_test, y_pred_3)
        ]
    }
    metric_eval_df = pd.DataFrame(metrics_eval)
    print(metric_eval_df)   

else : 
    mse = mean_squared_error(y_test, y_pred_1)
    mae = mean_absolute_error(y_test, y_pred_1)
    rmse = np.sqrt(mse)
    r2 = r2_score(y_test, y_pred_1)

    print(f"MSE : {mse:.4f}") 
    print(f"MAE : {mae:.4f}") 
    print(f"RMSE : {rmse:.4f}")
    print(f"R2 : {r2:.4f}") # part de variance expliqué par le model

print("\n")

# =============================================================================
# Importance des variables 
# =============================================================================

# Avec Score F = fréquence d’utilisation d’une variable comme critère de séparation dans les arbres de décision du model
if compare_model =="Y":
    xgb.plot_importance(xgb_model_1,title = "model 1")
    plt.show()
    xgb.plot_importance(xgb_model_2,title = "model 2")
    plt.show()
    xgb.plot_importance(xgb_model_3,title = "model 3")
    plt.show()
else :
    xgb.plot_importance(xgb_model_1,title = "model 1")
    plt.show()

# =============================================================================
# SHAP (SHapley Additive Explanations) 
# =============================================================================

if compare_model =="Y":
    explainer_1 = shap.Explainer(xgb_model_1, X_train)
    shap_values_1 = explainer_1(X_test)    
    explainer_2 = shap.Explainer(xgb_model_2, X_train)
    shap_values_2 = explainer_2(X_test)
    explainer_3 = shap.Explainer(xgb_model_3, X_train)
    shap_values_3 = explainer_3(X_test)

    shap.summary_plot(shap_values_1, X_test)
    shap.summary_plot(shap_values_2, X_test)
    shap.summary_plot(shap_values_3, X_test)
else : 
    explainer_1 = shap.Explainer(xgb_model_1, X_train)
    shap_values_1 = explainer_1(X_test)  
    shap.summary_plot(shap_values_1, X_test)

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
plt.scatter(X_train["prev_Mixo"], y_train, label="Valeurs train", color="grey", alpha=0.2)
plt.scatter(X_test["prev_Mixo"], y_test, label="Valeurs test réelles", color="blue", alpha=0.4)

# 1
plt.scatter(X_test["prev_Mixo"], y_pred_1, label="Valeurs prédites model 1", color="red", alpha=0.45)
plt.step(X_sorted, y_pred_sorted_1, where="post", color="red", linestyle="-", linewidth=1, alpha=0.5)

if compare_model =="Y":
    # 2
    plt.scatter(X_test["prev_Mixo"], y_pred_2, label="Valeurs prédites model 2", color="orange", alpha=0.45)
    plt.step(X_sorted, y_pred_sorted_2, where="post", color="orange", linestyle="-", linewidth=1, alpha=0.5)

    # 3
    plt.scatter(X_test["prev_Mixo"], y_pred_3, label="Valeurs prédites model 3", color="purple", alpha=0.45)
    plt.step(X_sorted, y_pred_sorted_3, where="post", color="purple", linestyle="-", linewidth=1, alpha=0.5)

plt.xlabel("Prévalence Mixotrophs")
plt.ylabel(f"Diversite : {choix_div}")
plt.legend()
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
    plt.subplot(4,2,2)
    plt.hist(residus_1, bins=30,color="grey")

    plt.subplot(4,2,3)
    plt.scatter(y_pred_2, residus_2, alpha=0.6)
    plt.axhline(y=0, color='r', linestyle='-', alpha=0.2)
    plt.subplot(4,2,4)
    plt.hist(residus_2, bins=30,color="grey")

    plt.subplot(4,2,5)
    plt.scatter(y_pred_3, residus_3, alpha=0.6)
    plt.axhline(y=0, color='r', linestyle='-', alpha=0.2)
    plt.subplot(4,2,6)
    plt.hist(residus_3, bins=30, color="grey")
    
    plt.subplot(4,2,7)
    X_train_const = sm.add_constant(X_train)  # Ajout d'une constante pour le modèle OLS
    model = sm.OLS(y_train, X_train_const).fit()
    influence = OLSInfluence(model)
    cooks_d = influence.cooks_distance[0]

    plt.stem(cooks_d, markerfmt=",")
    plt.ylabel("Cook’s Distance")

    plt.tight_layout()
    plt.show()

else:
    residus_1 = y_test - y_pred_1

    plt.scatter(y_pred_1, residus_1, alpha=0.6)
    plt.axhline(y=0, color='r', linestyle='-')
    plt.show()





























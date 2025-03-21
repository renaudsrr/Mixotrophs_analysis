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

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm

from scipy.stats import uniform,stats
from statsmodels.stats.outliers_influence import OLSInfluence

import xgboost as xgb
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

# Choix de paramétrage de modélisation
compare_model = input("(Vide = N) // Comparaison 3 model ? (Y/N)") or "N"

if compare_model =="Y":
    refaire_optimisation = input("(Vide = N) // Refaire l'optimisation hyperparams ? (Y/N) ") or "N"
    if refaire_optimisation =="Y":
        input_rd = input("(Vide = 100) // Nb essais optimisation hyperparam Randomised : ")
        n_trials_opt_Rd = int(input_rd) if input_rd.strip() != "" else 100

        input_op = input("(Vide = 80) // Nb essais optimisation hyperparam optuna : ")
        n_trials_opt_Op = int(input_op) if input_op.strip() != "" else 80

choix_div = input("(Vide = H) // Choix indice de diversite ? (rich_genus, H, 1-D, J) ") or "H"

suppr_ind_influents = input("(Vide = N) // Supprimer individus trop influents ? (Y/N) ")

# % données utilisées pour X_test et y_test
test_size = 0.1

print("\n")

# =============================================================================
# Data
# =============================================================================

data_python = pd.read_csv("Data/new_csv/data_python.csv",index_col=0)
VA_expl = data_python.drop(columns=["region","lake","type_prof","rich_genus","H","1-D","J"])

# Traitement des NaN (Provisioire) ============================================

VA_expl["biomass_tot_zoo"] = VA_expl["biomass_tot_zoo"].fillna(VA_expl["biomass_tot_zoo"].median())
VA_expl["sech.m"] = VA_expl["sech.m"].fillna(VA_expl["sech.m"].median())

VA_expl.to_csv("Data/new_csv/VA_expl_python.csv")

# =============================================================================
# Extraction data pour ML
# =============================================================================

y = data_python[choix_div] # Variable réponse
# Séparation avec jeu learning et test, ajout random_state = 42 si besoin pour cadrer résultats predi
X_train, X_test, y_train, y_test = train_test_split(VA_expl, y, test_size=test_size,random_state = 42) 

print("var train : ", np.var(y_train), " // var test : ", np.var(y_test))
print("\n")
# =============================================================================
# Supprimer individus trop influents
# =============================================================================

if suppr_ind_influents == "Y":
    X_train_const = sm.add_constant(X_train)
    model = sm.OLS(y_train, X_train_const).fit()

    # Calcul des distances de Cook
    influence = model.get_influence()
    cooks_d = influence.cooks_distance[0]

    n = len(X_train)
    threshold = 4 / n  # Seuil couramment utilisé

    influential_points = np.where(cooks_d > threshold)[0]
    influential_indices = X_train.index[influential_points]
    print("sample influents : ", list(influential_indices))
    X_train = X_train.drop(index=influential_indices)
    y_train = y_train.drop(index=influential_indices)

# =============================================================================
# Exploration data
# =============================================================================

plt.figure(figsize=(10, 8))

plt.scatter(X_train["prev_Mixo"],y_train,color="green",alpha=0.7,label="train")
plt.scatter(X_test["prev_Mixo"],y_test,color="green",alpha=0.3,label="test")

plt.legend()
plt.tight_layout()
plt.show()


# =============================================================================
# Optimisation hyperparamètres
# =============================================================================

if compare_model == "Y" and refaire_optimisation == "Y":

    print("opti Rd en cours...")

    # Randomise ==================================================================

    param_dist = {
        'n_estimators': [50, 100, 300, 500, 1000, 2000],
        # 'max_depth': [3, 6, 9, 12, 15, 20],
        # 'learning_rate': uniform(0.001, 0.3),  
        # 'subsample': uniform(0.5, 0.5),  
        'colsample_bytree': uniform(0.5, 0.5),
        # 'reg_lambda': uniform(0, 10),  
        # 'reg_alpha': uniform(0, 5),  
        'gamma' : uniform(0,0.5)
    }

    xgb_model_Rd = xgb.XGBRegressor(random_state=42)
    random_search = RandomizedSearchCV(xgb_model_Rd, param_dist, cv=5, scoring='r2', n_jobs=-1, n_iter=n_trials_opt_Rd)
    random_search.fit(X_train, y_train)

    # Bayesien Optuna ==================================================================

    print("opti Op en cours...")

    def objective(trial):
        params = {
            'n_estimators': trial.suggest_int('n_estimators', 50, 2000),
            # 'max_depth': trial.suggest_int('max_depth', 3, 10),
            # 'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.3, log=True),
            # 'subsample': trial.suggest_float('subsample', 0.5, 1.0),
            'colsample_bytree': trial.suggest_float('colsample_bytree', 0.5, 1.0),
            # 'reg_lambda': trial.suggest_float('reg_lambda', 0.0, 1.0),
            # 'reg_alpha': trial.suggest_float('reg_alpha', 0.0, 1.0),
            'gamma' : trial.suggest_float('gamma',0.0,0.5)
        }
        
        xgb_model_Op = xgb.XGBRegressor(**params)
        
        score = cross_val_score(xgb_model_Op, X_train, y_train, cv=5, scoring="r2",n_jobs=-1)
        return score.mean()

    study = optuna.create_study(direction='maximize')
    study.optimize(objective, n_trials=n_trials_opt_Op)

    hyperparams_Rd = pd.DataFrame.from_dict(random_search.best_params_, orient="index", columns=["Opti_Rd"])
    hyperparams_Op = pd.DataFrame.from_dict(study.best_params, orient="index", columns=["Opti_Op"])

    all_optimisation = pd.merge(hyperparams_Rd, hyperparams_Op, left_index=True, right_index=True, how="outer")
    ordre_hyperparams = ["learning_rate", "max_depth", "n_estimators", "colsample_bytree", "subsample", "reg_alpha", "reg_lambda", "gamma"]
    all_optimisation = all_optimisation.reindex(ordre_hyperparams)

    all_optimisation.to_csv(f"Data/new_csv/all_optimisation_{choix_div}.csv")

    print("\n")
    print(f"Optimisation hyperparams pour y = {choix_div}")
    print(all_optimisation)
    print(f"Meilleur R² Optuna : {study.best_value:.4f}")
    print("\n")

# =============================================================================
# Creation model
# =============================================================================

xgb_model_1 = xgb.XGBRegressor(
    random_state = 42,       # Graine aléatoire pour assurer la reproductibilité des résultats

    # grow_policy = "depthwise",   # Politique de croissance des arbres (par défaut : "depthwise", sinon "lossguide")

    # importance_type = "gain",   # Mesure l'importance des variables en fonction du gain apporté par chaque split

    # learning_rate = 0.1,    # Taux d'apprentissage (plus bas = apprentissage lent mais meilleure généralisation, "shrinkage")
    # max_depth = 15,           # Profondeur max des arbres (plus élevé = plus de complexité et risque d'overfitting)
    n_estimators = 500000,      # Nombre total d'arbres dans l'ensemble (plus élevé = modèle plus complexe et plus lent)

    colsample_bytree = 0.8,   # Proportion des colonnes (features) utilisées pour chaque arbre (ex: 80% des variables)
    # subsample = 0.8,          # Proportion des échantillons (lignes) utilisées pour chaque arbre (ex: 80% des lignes) 
    #                          # Permet d'ajouter de la variance et réduit le risque d'overfitting

    # reg_alpha = 0.0,         # Régularisation L1 (LASSO) : Encourage la réduction de certains coefficients à zéro (sparse model)
    # reg_lambda = 0.0,         # Régularisation L2 (Ridge) : Réduit la magnitude des coefficients pour éviter les valeurs extrêmes

    # min_child_weight = 4,    # Nombre minimal d'échantillons requis dans un nœud pour qu'il puisse être divisé 
                             # (plus haut = moins de splits, prévient l'overfitting)

    gamma = 0.2,            # Seuil minimum de réduction de la perte pour effectuer un split 
                             # (plus élevé = moins de splits inutiles)

    # objective = "reg:squarederror",   # Fonction de coût pour un problème de régression classique (MSE)

    # booster = "gbtree",      # Type de booster utilisé (gbtree = gradient boosting tree)

    eval_metric = "rmse"    # Indicateur d'évaluation utilisé (par défaut)
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

        # learning_rate = all_optimisation.loc["learning_rate", "Opti_Rd"],    
        # max_depth = int(all_optimisation.loc["max_depth", "Opti_Rd"]),          
        n_estimators = int(all_optimisation.loc["n_estimators", "Opti_Rd"]),    

        colsample_bytree = all_optimisation.loc["colsample_bytree", "Opti_Rd"],  
        # subsample = all_optimisation.loc["subsample", "Opti_Rd"],         

        # reg_alpha = all_optimisation.loc["reg_alpha", "Opti_Rd"],  
        # reg_lambda = all_optimisation.loc["reg_lambda", "Opti_Rd"],
        
        gamma = all_optimisation.loc["gamma", "Opti_Rd"],  
        )

    xgb_model_3 = xgb.XGBRegressor(
        random_state = 42,        
        # grow_policy = "depthwise",
        # importance_type = "gain",
        
        # learning_rate = all_optimisation.loc["learning_rate", "Opti_Op"],    
        # max_depth =  int(all_optimisation.loc["max_depth", "Opti_Op"]),          
        n_estimators =  int(all_optimisation.loc["n_estimators", "Opti_Op"]),    

        colsample_bytree =  all_optimisation.loc["colsample_bytree", "Opti_Op"],  
        # subsample =  all_optimisation.loc["subsample", "Opti_Op"],         

        # reg_alpha =  all_optimisation.loc["reg_alpha", "Opti_Op"],  
        # reg_lambda =  all_optimisation.loc["reg_lambda", "Opti_Op"],
        
        gamma = all_optimisation.loc["gamma", "Opti_Op"],  
        )
                                             
# =============================================================================
# Learning model
# =============================================================================

if compare_model=="Y":
    xgb_model_1.fit(X_train,y_train)
    xgb_model_2.fit(X_train,y_train)
    xgb_model_3.fit(X_train,y_train)
else:
    xgb_model_1.fit(X_train,y_train, verbose=False)

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
    xgb.plot_importance(xgb_model_2,title = "model 2")
    xgb.plot_importance(xgb_model_3,title = "model 3")
    plt.show()
else :
    xgb.plot_importance(xgb_model_1,title = "model 1")
    plt.show()

fi = pd.DataFrame(data=xgb_model_1.feature_importances_,
    index=xgb_model_1.feature_names_in_, columns=["importance"]).sort_values("importance")

# =============================================================================
# SHAP (SHapley Additive Explanations) 
# =============================================================================

# Avec valeur SHAP = indique dans quelle mesure une variable influence la prédiction du modèle. (Shap positif => influence fit positif). couleur = valeut de la variable

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
    X_train_const = sm.add_constant(X_train) 
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

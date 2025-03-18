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

import xgboost as xgb
import shap

from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.inspection import PartialDependenceDisplay

# =============================================================================
# Data
# =============================================================================

data_phyto = pd.read_csv("Data/new_csv/data_phyto.csv") 
data_zoo = pd.read_csv("Data/new_csv/data_zoo.csv")

# DATA AVEC : 1ère VA Reponse = 4 métriques de diversité phyto & 1 VA explicative = prev_Mixo
summary_data_phyto = pd.read_csv("Data/new_csv/summary_data_phyto.csv") 
#  DATA AVEC : 2ème VA Reponse = 4 métriques de diversité zoo
summary_data_zoo = pd.read_csv("Data/new_csv/summary_data_zoo.csv") 

# VA explicatives environementales 
env = pd.read_csv("Data/new_csv/env_modif.csv") 

# pour strat mixotrophs
info_genus = pd.read_csv("Data/new_csv/info_genus_zoo.csv") 

# fusionner prev_Mixo avec Env pour data_explicatifs
data_explicatif = pd.concat([env,summary_data_phyto["prev_Mixo"]],axis=1)
# print(data_explicatif.dtypes)
data_explicatif = data_explicatif.drop(columns=["Sample","region","lake","type_prof","doy","yr"])


# =============================================================================
# Extraction data pour ML
# =============================================================================

choix_div = "1-D"
y = summary_data_phyto[choix_div] # Variable réponse

# Séparation avec jeu learning et test
X_train, X_test, y_train, y_test = train_test_split(data_explicatif, y, test_size=0.2, random_state=42)

# =============================================================================
# Creation model
# =============================================================================

xgb_model = xgb.XGBRegressor(
    n_estimators = 100,     # nb arbres
    learning_rate = 0.01,    # Plus il est bas, plus le modèle apprend lentement mais en généralisant mieux => Shrinkage
    max_depth = 10,           # Profondeur max des arbres => traux haut = overfitting
    subsample = 0.8,         # Prend 80 % du jeu learning pour chaque arbrew => evite overfitting
    colsample_bytree = 0.8,  # Prend 80 % des variables dans chaque arbre
    random_state = 42        # graine aleatoire pour reproductibilite
    )

# =============================================================================
# Learning model
# =============================================================================

xgb_model.fit(X_train,y_train)

# =============================================================================
# Prediction
# =============================================================================

y_pred = xgb_model.predict(X_test)

# =============================================================================
# Evaluation
# =============================================================================

# Avec MSE
mse = mean_squared_error(y_test, y_pred)

# Avec R2
r2 = r2_score(y_test, y_pred)

print(f"MSE : {mse:.4f}") # précision des prédictions y <-> y^
print(f"R2 : {r2:.4f}") # part de variance expliqué par le model

# =============================================================================
# Importance des variables 
# =============================================================================

xgb.plot_importance(xgb_model)
plt.show()

# =============================================================================
# SHAP (SHapley Additive Explanations) 
# =============================================================================

explainer = shap.Explainer(xgb_model, X_train)
shap_values = explainer(X_test)

shap.summary_plot(shap_values, X_test)

# =============================================================================
# Relation X-y
# =============================================================================

sorted_indices = np.argsort(X_test["prev_Mixo"])
X_sorted = X_test["prev_Mixo"].iloc[sorted_indices]
y_pred_sorted = y_pred[sorted_indices]

plt.figure(figsize=(8,6))
plt.scatter(X_test["prev_Mixo"], y_test, label="Valeurs réelles", color="blue", alpha=0.7)
plt.scatter(X_test["prev_Mixo"], y_pred, label="Valeurs prédites", color="red", alpha=0.7)
plt.plot(X_sorted, y_pred_sorted, color="red", linestyle="-", linewidth=1)
plt.xlabel("Prévalence Mixotrophs")
plt.ylabel(f"Diversite : {choix_div}")
plt.legend()
plt.show()

plt.figure(figsize=(8,6))
plt.scatter(X_test["prev_Mixo"], X_test["ctch.m2"], label="Valeurs réelles", color="blue", alpha=0.7)
plt.scatter(X_test["prev_Mixo"], X_test["ctch.m2"], label="Valeurs prédites", color="red", alpha=0.7)
plt.xlabel("Prévalence Mixotrophs")
plt.ylabel(f"Diversite : {choix_div}")
plt.legend()
plt.show()



# =============================================================================
# Optimisation hyperparamètres
# =============================================================================

# param_grid = {
#     'n_estimators': [50,100, 300, 500,1000],
#     'max_depth': [3, 6, 9, 15],
#     'learning_rate': [0.005,0.01, 0.05, 0.1],
#     'subsample': [0.7, 0.8, 0.9],
#     'colsample_bytree': [0.7, 0.8, 0.9]
# }

# xgb_model = xgb.XGBRegressor(random_state=42)
# grid_search = GridSearchCV(xgb_model, param_grid, cv=5, scoring='r2', n_jobs=-1,verbose=3)

# grid_search.fit(X_train, y_train)
# print("Meilleurs hyperparamètres :", grid_search.best_params_)



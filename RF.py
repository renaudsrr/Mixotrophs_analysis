"""
Script building a model [...]

Author : SERRE Renaud
Creation 12.03.2025

History of modification
12.03.2025 : Creation of the script. Import data 
13.03.2025 : Some test of RF

"""

# =============================================================================
# Packages
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import xgboost as xgb
import shap

# =============================================================================
# Data
# =============================================================================

data_phyto = pd.read_csv("Data/new_csv/data_phyto.csv") 
data_zoo = pd.read_csv("Data/new_csv/data_zoo.csv")

summary_data_phyto = pd.read_csv("Data/new_csv/summary_data_phyto.csv")
summary_data_zoo = pd.read_csv("Data/new_csv/summary_data_zoo.csv")

info_genus = pd.read_csv("Data/new_csv/info_genus_zoo.csv")

env_biol = pd.read_csv("Data/EnvData/Biol.csv",sep=";")
env_CM = pd.read_csv("Data/EnvData/Carbon_Metabol.csv",sep=";")
env_chim = pd.read_csv("Data/EnvData/Chim.csv",sep=";")
env_GPP = pd.read_csv("Data/EnvData/GPP.csv",sep=";")
env_phy = pd.read_csv("Data/EnvData/Phy.csv",sep=";")

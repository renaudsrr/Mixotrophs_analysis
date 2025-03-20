# **TITRE**

## **Scripts by SERRE Renaud**  
**Internship UQAM-GRIL under the supervision of BEISNER Beatrix**

This project studies the link between mixotrophic prevalence and phytoplankton diversity in lakes. The work includes data processing in R and machine learning modeling in Python.

---

## **1. Data Processing (`exploration_data_CARBBAS.R`)**

This R script cleans and prepares data on biodiversity and environmental factors.

### **Main Steps**  
- **Data Import & Cleaning:** Loads and organizes data on phytoplankton, zooplankton, and the environment.
- **Biological Metrics Calculation:**  
  - Computes total biomass for phytoplankton and zooplankton.  
  - Measures species richness & calculates diversity indices (Shannon, Simpson, and Pielou’s evenness).  
  - Estimates mixotrophic species presence and biomass.  
- **Environmental Data Integration:** Merges biological and environmental data.
- **Data Export:** Saves processed data in the `new_csv/` folder for machine learning.

---

## **2. Machine Learning (`ML.py`)**

This Python script builds models to study relation between prev_Mixo and parameter of diversity taking into account environement parameters.

### **Main Steps**  
- **Exploratory Analysis:** Plots relationships between mixotrophic prevalence and diversity.
- **Feature Processing & Normalization:** Handles missing values and scales data if needed.
- **Model Selection & Optimization:**  
  - Uses **XGBoost** for regression.  
  - Tunes hyperparameters with **RandomizedSearchCV** and **Optuna**.
- **Model Evaluation:**  
  - Calculates metrics of evaluation scores.  
  - Creates **SHAP** plots to show feature importance.  
  - Analyzes model errors.
- **Predictions & Visualizations:**  
  - Compares real vs. predicted values.  
  - Evaluates model performance.

---
This work helps understand how mixotrophic prevalence affect phytoplankton diversity in lakes.

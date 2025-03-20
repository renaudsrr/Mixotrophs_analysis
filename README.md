# **Analysis of Biodiversity and Environmental Factors in Lakes**

## **Scripts by SERRE Renaud**  
**Internship UQAM-GRIL under the supervision of BEISNER Beatrix**

This project explores the relationship between mixotrophic prevalence and phytoplankton diversity using ecological data from multiple lakes. The workflow consists of data processing in R and machine learning modeling in Python.

---

## **1. Data Processing (`exploration_data_CARBBAS.R`)**

This R script performs data extraction, cleaning, and computation of key biodiversity and environmental metrics.

### **Main Features**  
- **Data Import & Cleaning:** Loads raw datasets on phytoplankton, zooplankton, and environmental parameters, ensuring consistency across sources.
- **Biological Metrics Calculation:**  
  - Computes total biomass for phytoplankton and zooplankton.  
  - Estimates species richness and taxonomic classification.  
  - Calculates diversity indices (Shannon, Simpson, Pielou’s evenness).  
  - Assesses mixotrophic species abundance and biomass.  
- **Environmental Data Integration:** Merges biological and environmental datasets for further analysis.
- **Data Export:** Processed data is saved in the `new_csv/` folder for use in machine learning.

---

## **2. Machine Learning Modeling (`ML.py`)**

This Python script builds predictive models to analyze the relationship between environmental variables and biodiversity metrics.

### **Key Components**  
- **Data Preparation:** Loads and preprocesses cleaned data from `new_csv/`.
- **Exploratory Analysis:** Visualizes relationships between mixotrophic prevalence and diversity indices.
- **Feature Engineering & Normalization:** Handles missing values and applies standardization when needed.
- **Model Selection & Optimization:**  
  - Implements **XGBoost** for regression analysis.  
  - Compares multiple models with hyperparameter optimization using **RandomizedSearchCV** and **Optuna**.
- **Model Evaluation:**  
  - Computes **MSE, RMSE, MAE, and R²** scores.  
  - Generates **SHAP** plots for feature importance analysis.  
  - Examines residuals to assess model quality.
- **Predictions & Visualizations:**  
  - Plots relationships between observed and predicted values.  
  - Compares performance of optimized models.

---

## **Workflow Summary**
1. **Data Cleaning & Metrics Calculation** (R) → Export to `new_csv/`
2. **Exploratory Data Analysis** (Python)
3. **Machine Learning Modeling & Optimization** (Python)
4. **Model Interpretation & Evaluation**

This structured approach allows for a robust assessment of the factors influencing phytoplankton diversity and mixotrophic prevalence in lake ecosystems.

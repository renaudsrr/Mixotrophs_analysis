# Aquatic Biodiversity Modeling Project

This project focuses on harmonizing, exploring, and modeling planktonic biodiversity in North American lakes using ecological statistics and machine learning. It combines R scripts for ecological data processing and analysis with Python scripts for predictive modeling and user interaction.

---

## General Project Workflow

1. **Raw Data Processing**  
   - **Scripts**: `1_Processing_CARBBAS.R`, `2_Processing_LP_NLA.R`, `3_Processing_IISD_ELA.R`, `Processing_bact_LP.R`  
   - These scripts clean, aggregate, and format raw phytoplankton, zooplankton, and bacteria datasets from different lake monitoring programs (CARBBAS, LP-NLA, ELA).

2. **Data Harmonization**  
   - **Script**: `Harmonisation_data.R`  
   - Merges and harmonizes the cleaned datasets into a unified format suitable for cross-study comparison and downstream analyses.

3. **Diversity and Taxonomic Analyses**  
   - **Scripts**: `Div_BETA.R`, `Tanonomic_analysis_all_data.R`, `infos_taxonomy.R`  
   - Computes alpha and beta diversity metrics, summarizes taxonomic distributions, and links genus/species names to higher-level taxonomy.

4. **Statistical Modeling**  
   - **Scripts**: `GAM_model.R`, `Verif_PB_circularite.R`  
   - Fits Generalized Additive Models (GAMs) to investigate relationships between diversity metrics and environmental drivers.
   - Includes verification of circularity or bias issues in modeling.

5. **Spatial Analysis**  
   - **Script**: `Spatial_analysis.R`  
   - Explores spatial patterns in biodiversity across lakes, using geospatial coordinates to detect spatial structure or clustering.

6. **Machine Learning Interface**  
   - **Scripts**: `ML_v13.py`, `open_interface.py`  
   - Python-based tools for machine learning modeling, integrating algorithms like XGBoost and Random Forest.
   - `open_interface.py` provides a Tkinter graphical interface for configuring datasets, models, and evaluation options.

---

## Integration and Purpose

The pipeline is designed to:
- Prepare high-quality biodiversity datasets across different monitoring programs.
- Analyze ecological diversity through both classical and spatial statistics.
- Use advanced regression and machine learning techniques to predict diversity metrics.
- Provide interactive tools for reproducible and flexible analysis workflows.

---

## Technologies Used

- **R**: ecological data wrangling, statistical modeling (vegan, mgcv), and plotting.
- **Python**: machine learning (XGBoost, scikit-learn), interface development (Tkinter).
- **GeoJSON** and spatial tools for mapping lake locations.

---

## File Overview

| File Name                     | Purpose                                               |
|------------------------------|--------------------------------------------------------|
| 1_Processing_CARBBAS.R       | Prepares CARBBAS dataset                              |
| 2_Processing_LP_NLA.R        | Prepares LakePulse + NLA dataset                      |
| 3_Processing_IISD_ELA.R      | Prepares ELA dataset                                  |
| Processing_bact_LP.R         | Processes bacteria data from LP                       |
| Harmonisation_data.R         | Harmonizes all datasets into a common format          |
| Div_BETA.R                   | Computes beta diversity and plots                     |
| infos_taxonomy.R             | Adds taxonomic metadata                               |
| Tanonomic_analysis_all_data.R| Summary of taxonomic structure across datasets        |
| GAM_model.R                  | Fits GAMs to model relationships with diversity        |
| Verif_PB_circularite.R       | Checks for modeling circularity issues                |
| Spatial_analysis.R           | Performs geospatial biodiversity analysis             |
| ML_v13.py                    | Runs machine learning models                          |
| open_interface.py            | GUI for configuring and launching ML models           |

---

## Author

Renaud Serre – Master's intern at UQAM-GRIL  
Supervised by Beatrix Beisner

---

# TITRE

## Script by SERRE Renaud

**Internship UQAM-GRIL under the supervision of BEISNER Beatrix**

### Data analysis and manipulation => à modifier à la fin 

The R script (`exploration_data.R`) processes and analyzes biological and environmental data from different lakes. It extracts key metrics related to phytoplankton, zooplankton, and environmental parameters. The processed data is later used in an XGBoost model in the `RF.py` script.

### Main Features

#### 1. **Summary Statistics Calculation**

- Computes total biomass for phytoplankton and zooplankton.
- Calculates species richness for each lake.
- Merges taxonomic and functional information for genus classification strategy.

#### 2. **Biodiversity Metrics**

- Calculates Shannon, Simpson diversity indices and Pielou’s evenness
- Estimates mixotrophic species abundance and biomass.

#### 3. **Export Processed Data**

- Saves cleaned and summarized datasets to `new_csv/` folder :

### RF

- The processed data will be used for machine learning analysis with XGBoost in `RF.py`.

---

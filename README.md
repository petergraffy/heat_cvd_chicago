# Modeling Excess Cardiovascular Mortality During Extreme Heat in Chicago (2010–2022)

This repository contains all analytic code, data processing scripts, and derived datasets used in our study examining how extreme heat affects cardiovascular disease (CVD) mortality across Chicago community areas from 2010–2022. The project integrates mortality records, high-resolution temperature and humidity data, ACS sociodemographic estimates, and several environmental exposures (PM₂.₅, NO₂, NDVI) to generate multi-scale mortality models and community-specific heat thresholds.

## 📁 Repository Structure

### Data (CSV / GeoJSON)
These files are cleaned, analysis-ready datasets used to build annual, monthly, weekly, and daily modeling panels.

| File | Description |
|------|-------------|
| **acs_ca_2010–2022.csv** | Community-area ACS estimates (demographics, socioeconomic indicators). |
| **chicago_acs_2010–2022.csv** | Additional ACS variables harmonized for citywide analyses. |
| **CA_temps_rh_90–23.csv** | Daily temperature & humidity (Daymet + NOAA) aggregated to Chicago community areas. |
| **areas_tracts2010_joined.csv** | Census tract ↔ community area spatial crosswalk (2010 TIGER). |
| **comm_areas.geojson** | Chicago community area shapefile used for all spatial joins & mapping. |
| **no2_panel_community_year.csv** | Annual community-area NO₂ concentrations from multi-source datasets. |
| **pm25_panel_community_year.csv** | Annual community-area PM₂.₅ concentrations (1998–2023 coverage). |
| **no2_county_year.csv / pm25_county_year.csv** | CONUS county-year pollutant datasets for cross-validation and sensitivity analyses. |
| **ndvi_year_ca.csv** | Community-area NDVI (vegetation greenness) time-series. |

## 🧰 Code

### Main Analysis Pipeline
- **run_heat_cvd_all_scales.r**  
  Master script that builds daily, monthly, and annual analytic panels and runs all modeling steps:
  - GAMs and quasi-Poisson models  
  - Distributed-lag sensitivity analyses  
  - Community-specific threshold detection  
  - Excess mortality estimation at multiple temporal scales  
  - Figure generation (publication-ready outputs)

### Environmental Data Acquisition
- **ndvi_download.r / ndvi_download_agg.r**  
  Scripts to download, mosaic, aggregate, and extract NDVI for each community area.

## 🔬 Analytic Overview

This project constructs a high-resolution, longitudinal dataset integrating:

### 1. Mortality Data
- Chicago Department of Public Health death records, 2010–2022  
- CVD outcomes identified using ICD-10 I-codes  
- Daily counts aggregated at the community-area level

### 2. Climate Data
- Daymet 1-km max/min/mean temperature  
- NOAA relative humidity  
- Temporal scales: **daily**, **daily moving averages**, **monthly**, **warm-season**, and **annual**

### 3. Environmental Exposures
- PM₂.₅ (1998–2023) and NO₂ (2005–2020)  
- NDVI (annual greenness)  
- Spatial resolution harmonized to community areas

### 4. Demographics
- 5-year ACS indicators: population, age, race/ethnicity, education, income, insurance, employment, etc.

### 5. Modeling Framework
- GAMs with smooth terms for temperature, humidity, age, income, and time  
- Quasi-Poisson models for rates per 100,000  
- Community-specific heat thresholds using nonlinear exposure-response curves  
- Excess mortality estimation under observed vs counterfactual (non-heat) scenarios  
- Sensitivity analyses excluding COVID years  

## 📊 Outputs

Running the pipeline produces:

- **Heat–CVD exposure–response curves** (daily, monthly, annual)  
- **Community-specific heat thresholds**  
- **Excess deaths during extreme heat events**  
- **Spatial maps** of thresholds, excess deaths, and annual trends  
- Sensitivity analyses (e.g., excluding COVID, alternative lag structures)  
- Cleaned and reproducible analytic panels  

## 🚀 Reproducibility

To run the entire analysis:

```r
source("run_heat_cvd_all_scales.r")
```

Required R packages include:  
`tidyverse`, `sf`, `mgcv`, `data.table`, `tmap`, `lubridate`, `vroom`, and others listed at the top of the master script.

## 🔍 Citation

If you use this repository, please cite:

**Graffy P., et al. Modeling Excess Cardiovascular Mortality During Extreme Heat Across Chicago Community Areas, 2010–2022.**  
*Manuscript in preparation.*

## 📫 Contact

For questions or collaboration requests:

**Peter Graffy, PhD, MPH**  
University of Chicago  
🌐 https://petergraffy.github.io/

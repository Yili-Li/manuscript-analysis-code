# Analysis code: Association between the weight-adjusted-waist index (WWI) and breast cancer prevalence in U.S. adult women (NHANES 1999-2016)

Core analysis code examining the association between the weight-adjusted-waist
index (WWI) and breast cancer prevalence in U.S. adult women, based on NHANES
1999-2016 public-use data.

## Contents

| Folder | Script | Analysis |
|---|---|---|
| 01_weighted_logistic | weighted_logistic_models.R | Survey-weighted binary logistic regression: WWI continuous and in quartiles, three models with staged covariate adjustment |
| 02_rcs | rcs_analysis.R | Restricted cubic splines (3 knots) with adjusted model; exploratory inflection point with bootstrap CI |
| 03_subgroup_interaction | subgroup_forest_interaction.R | Prespecified subgroup analyses (crude within-stratum ORs) and unadjusted WWI x subgroup interaction tests with Benjamini-Hochberg FDR |
| 04_psm_sensitivity | psm_matching_balance.R | 1:1 propensity-score matching (nearest neighbor, without replacement, caliper 0.2 SD of the logit PS); post-match minimal-adjustment models |
| 04_psm_sensitivity | love_plot_smd.R | Covariate balance before/after matching: variable-level standardized mean differences (Love plot) |
| 05_xgboost_shap | xgboost_shap_pipeline.R | XGBoost variable importance and SHAP values on the full analytic cohort; ten 2,000-participant subsample replicates |
| 05_xgboost_shap | xgboost_cv_nrounds.R | Post-hoc 10-fold cross-validation: AUC vs number of boosting iterations |
| 06_hormones | hormone_weighted_models.R | Survey-weighted linear regression of sex-hormone biomarkers on WWI (women >= 20 years, NHANES 2013-2016) |

## Software

R 4.5.1. Packages: survey, rms, ggrcs, jstable, Matching, tableone, xgboost,
shapviz, pROC, ggplot2.

## Data

Public-use NHANES data (1999-2016): https://www.cdc.gov/nchs/nhanes/

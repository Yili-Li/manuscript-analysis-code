# ============================================================================
# XGBoost variable importance and SHAP analysis
# (manuscript Section 2.10; Figure 2; Table S6)
#
# Descriptive variable-importance ranking on the full analytic cohort
# (n = 6,782; 208 breast-cancer cases). Predictors: age (dichotomized at
# 60 years), race, PIR, marital status, alcohol use, hypertension, diabetes,
# CHD, heart failure, and WWI (entered continuously). The default regression
# objective (reg:squarederror) is used, with a learning rate of 0.05 and
# 6 boosting iterations; explicit train/test splitting and hyperparameter
# tuning are not performed because the purpose is descriptive importance
# extraction in a rare-outcome sample rather than prediction. Survey weights
# are not applied. SHAP values are computed for the full cohort via
# predcontrib (shapviz) and interpreted on the model raw output scale. As a
# robustness check, the importance ordering is re-estimated in ten
# independent 2,000-participant subsamples, both with the fixed model on
# re-drawn explanation sets and with models refitted on the subsamples.
# ============================================================================

suppressPackageStartupMessages({ library(xgboost); library(shapviz); library(patchwork) })

## ---- Input: derived analytic cohort ------------------------------------------
## X: feature data frame (n = 6,782 x 10) coded as follows:
##   Age (0/1 = < 60 / >= 60), Race (1-5), PIR (1-3), Marital (1/0),
##   Alcohol (1/0), Hypertension (1/0), Diabetes (1/0), CHD (1/0),
##   HeartFailure (1/0), WWI (continuous);
## y: outcome vector (1 = breast cancer, 0 = control).
# load("analytic_cohort_xgboost.RData")   # <- placeholder

xvars <- colnames(X)

## ---- Train (learning rate 0.05, 6 boosting iterations, default objective) ----
set.seed(20241115)
fit <- xgb.train(params = list(learning_rate = 0.05),
                 data   = xgb.DMatrix(data.matrix(X), label = y),
                 nrounds = 6)

## ---- Full-cohort SHAP ---------------------------------------------------------
shp_full <- shapviz(fit, X_pred = data.matrix(X), X = X)   # predcontrib
imp_full <- sort(colMeans(abs(shp_full$S[, xvars, drop = FALSE])),
                 decreasing = TRUE)
print(round(imp_full, 5))          # cf. Section 3.3 / Figure 2 panel A

## ---- Robustness: ten independent 2,000-participant subsamples ----------------
rank_explain <- rank_refit <- matrix(NA_real_, length(xvars), 10,
                                     dimnames = list(xvars, paste0("s", 1:10)))
for (i in 1:10) {
  set.seed(20241115 + i)
  idx <- sample(nrow(X), 2000)
  Xe  <- X[idx, ]

  ## (a) fixed model, re-drawn explanation set
  si <- shapviz(fit, X_pred = data.matrix(Xe), X = Xe)
  rank_explain[, i] <- rank(-colMeans(abs(si$S[, xvars, drop = FALSE])))

  ## (b) model refitted on the subsample, explained on the full cohort
  mf <- xgb.train(params = list(learning_rate = 0.05),
                  data   = xgb.DMatrix(data.matrix(Xe), label = y[idx]),
                  nrounds = 6)
  sf <- shapviz(mf, X_pred = data.matrix(X), X = X)
  rank_refit[, i] <- rank(-colMeans(abs(sf$S[, xvars, drop = FALSE])))
}
stab <- data.frame(Feature = xvars,
                   rank_full          = rank(-imp_full)[xvars],
                   explain_rank_min   = apply(rank_explain, 1, min),
                   explain_rank_max   = apply(rank_explain, 1, max),
                   refit_rank_min     = apply(rank_refit, 1, min),
                   refit_rank_max     = apply(rank_refit, 1, max))
stab[order(stab$rank_full), ]      # cf. Supplementary Table S6

## ---- Figure 2: importance (with numbers) + beeswarm summary -------------------
p_imp <- sv_importance(shp_full, show_numbers = TRUE)
p_bee <- sv_importance(shp_full, kind = "beeswarm")
fig <- p_imp + p_bee
ggsave("Figure2_xgboost_importance_shap.pdf", fig, width = 10, height = 6, units = "in")
ggsave("Figure2_xgboost_importance_shap.png", fig, width = 10, height = 6,
       units = "in", dpi = 300)

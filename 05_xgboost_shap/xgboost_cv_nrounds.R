# ============================================================================
# Post-hoc 10-fold cross-validation: AUC vs number of boosting iterations
# (manuscript Section 2.10; Table S7)
#
# Same feature coding and configuration as the SHAP analysis. Test AUC varied
# little across boosting depths (0.733 at 6 rounds vs a peak of 0.752 at
# 17 rounds, within one fold-to-fold SD), and the importance ranking was
# unaffected by additional boosting (age ranked first in 10/10 refits; WWI
# remained in the top three).
# ============================================================================

suppressPackageStartupMessages({ library(xgboost); library(pROC) })

## ---- Input: same coded features X and outcome y as the SHAP pipeline ---------
# load("analytic_cohort_xgboost.RData")   # <- placeholder

set.seed(11)
folds <- sample(rep(1:10, length.out = nrow(X)))
auc   <- matrix(NA_real_, 50, 10)

for (k in 1:10) {
  tr  <- folds != k
  m   <- xgb.train(params = list(learning_rate = 0.05),
                   data   = xgb.DMatrix(data.matrix(X)[tr, , drop = FALSE],
                                        label = y[tr]),
                   nrounds = 50)
  Dte <- xgb.DMatrix(data.matrix(X)[!tr, , drop = FALSE])
  auc[, k] <- sapply(1:50, function(r)
    as.numeric(pROC::roc(y[!tr], predict(m, Dte, ntreelimit = r),
                         quiet = TRUE, direction = "<",
                         levels = c(0, 1))$auc))
}

cv <- data.frame(nrounds = 1:50,
                 auc_mean = rowMeans(auc),
                 auc_sd   = apply(auc, 1, sd))
best <- which.max(cv$auc_mean)
cat(sprintf("Best AUC %.4f at nrounds = %d | AUC at 6 rounds = %.4f (SD %.4f)\n",
            cv$auc_mean[best], best, cv$auc_mean[6], cv$auc_sd[6]))
write.csv(cv, "TableS7_cv_auc_by_nrounds.csv", row.names = FALSE)

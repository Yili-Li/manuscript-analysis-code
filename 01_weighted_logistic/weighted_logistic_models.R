# ============================================================================
# Survey-weighted logistic regression: WWI and breast cancer prevalence
# (manuscript Section 2.6; Table 2)
#
# Design: NHANES 1999-2016, eight pooled cycles; MEC examination weights
#         (WTMEC2YR divided by the number of pooled cycles), masked strata
#         (SDMVSTRA), primary sampling units (SDMVPSU), Taylor-series
#         linearization.
# Models: Model 1 unadjusted;
#         Model 2 adjusted for age, race, and poverty-to-income ratio (PIR);
#         Model 3 further adjusted for hypertension, diabetes, and CHD
#         (i.e., Model 3 includes age, race, PIR, hypertension, diabetes, CHD).
#         Age is used as a continuous variable in adjusted models.
#         WWI enters as a continuous variable and in quartiles (Q1 reference).
# Reporting: odds ratios with profile-likelihood 95% CIs (rounded to three
#         decimals); P-values are Wald tests.
# ============================================================================

suppressMessages(library(survey))

## ---- Input: derived analytic cohort (n = 6,782; 208 breast-cancer cases) ----
## Required columns: breast (0/1), Age (years, continuous), Race (5 categories),
## PIR (low / middle / high), Hypertension, Diabetes, CHD (0/1),
## WWI (cm/sqrt(kg)), and design variables WTMEC2YR, SDMVSTRA, SDMVPSU.
# load("analytic_cohort_main.RData")   # <- derived analytic file (placeholder)

cohort$WWI_Q <- cut(cohort$WWI,
                    breaks = quantile(cohort$WWI),
                    include.lowest = TRUE,
                    labels = c("Q1", "Q2", "Q3", "Q4"))   # Q1 = reference

## ---- Survey design ----------------------------------------------------------
des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA,
                 weights = ~WTMEC2YR/8, nest = TRUE, data = cohort)

## ---- Model formulas (continuous and quartile WWI) ---------------------------
cov_m2 <- "Age + Race + PIR"
cov_m3 <- "Age + Race + PIR + Hypertension + Diabetes + CHD"
forms  <- list(
  M1_cont = "breast ~ WWI",
  M2_cont = paste("breast ~ WWI +", cov_m2),
  M3_cont = paste("breast ~ WWI +", cov_m3),
  M1_Q    = "breast ~ WWI_Q",
  M2_Q    = paste("breast ~ WWI_Q +", cov_m2),
  M3_Q    = paste("breast ~ WWI_Q +", cov_m3))

fits <- lapply(forms, function(f)
  svyglm(as.formula(f), design = des, family = quasibinomial()))

## ---- Extract WWI effects: profile-likelihood 95% CI, Wald P -----------------
summarize_term <- function(fit, term) {
  est <- coef(fit)[term]
  ci  <- confint(fit, term, method = "likelihood")     # profile likelihood
  p   <- summary(fit)$coef[term, "Pr(>|t|)"]           # Wald test
  data.frame(OR = exp(est), lower = exp(ci[1]), upper = exp(ci[2]), P = p)
}
lapply(fits[c("M1_cont", "M2_cont", "M3_cont")], summarize_term, term = "WWI")
lapply(fits[c("M1_Q", "M2_Q", "M3_Q")], summarize_term, term = "WWI_QQ4")  # Q4 vs Q1

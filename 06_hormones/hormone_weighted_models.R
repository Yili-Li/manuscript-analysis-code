# ============================================================================
# Survey-weighted linear regression: WWI and sex-hormone biomarkers
# (manuscript Section 2.2; Table 3)
#
# Cohort: women aged >= 20 years, NHANES 2013-2014 and 2015-2016 cycles
# (serum estradiol, total testosterone, and SHBG were measured by
# isotope-dilution LC-MS/MS only in these cycles). After merging hormone
# laboratory data with the main analytic cohort (intermediate n = 4,006) and
# excluding biologically implausible concentrations among women
# (estradiol > 1,000 pg/mL or total testosterone > 200 ng/dL; n = 32, all
# without breast cancer), the final hormone cohort comprises n = 3,974
# (117 breast-cancer cases and 3,857 participants without breast cancer).
# Model: Outcome ~ WWI + age + race + PIR + hypertension + diabetes + CHD
# (per 1-unit increase in WWI, original concentration scale).
# Weights: MEC examination weights (WTMEC2YR divided by 2 for two pooled
# cycles) with Taylor-series linearization (SDMVSTRA / SDMVPSU, nest = TRUE).
# CIs use the design degrees of freedom (t distribution); P-values are Wald
# tests. The Control column is restricted to the 3,857 cancer-free women.
# ============================================================================

suppressMessages(library(survey))

## ---- Input: hormone analysis cohort (intermediate cohort, n = 4,006) --------
## breast (Breast cancer / Control), WWI, age, race, pir, htn, dm, chd,
## E2 (pg/mL), TT (ng/dL), SHBG (nmol/L), and design variables
## WTMEC2YR, SDMVSTRA, SDMVPSU.
# load("analytic_cohort_hormones.RData")   # <- placeholder

d <- d[!with(d, E2 > 1000 | TT > 200), ]   # biologically implausible values

## ---- Survey design (Overall + cancer-free subset) ----------------------------
des      <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA,
                      weights = ~I(WTMEC2YR/2), nest = TRUE, data = d)
des_ctrl <- subset(des, breast == "Control")

## ---- One weighted regression: beta, SE, 95% CI (design df), Wald P, n --------
run_one <- function(design, outcome) {
  m  <- svyglm(as.formula(paste(outcome,
                                "~ WWI + age + race + pir + htn + dm + chd")),
               design = design)
  s  <- summary(m)$coef
  df <- degf(design)
  ci <- s["WWI", 1] + c(-1, 1) * qt(0.975, df) * s["WWI", 2]
  data.frame(beta = s["WWI", 1], SE = s["WWI", 2],
             lower = ci[1], upper = ci[2],
             P = 2 * pt(-abs(s["WWI", 1] / s["WWI", 2]), df),
             n = nrow(design$variables))
}

outcomes <- c(E2 = "E2", TT = "TT", SHBG = "SHBG",
              "TT/E2" = "I(TT/E2)", "TT/SHBG" = "I(TT/SHBG)")

res_overall <- lapply(outcomes, run_one, design = des)        # Overall (n = 3,974)
res_control <- lapply(outcomes, run_one, design = des_ctrl)   # Control (n = 3,857)
res_overall; res_control

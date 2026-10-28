# ============================================================================
# Propensity-score matching and post-match models
# (manuscript Sections 2.5 and 2.9; Tables S1 and S5)
#
# Propensity scores are estimated with a logistic regression model including
# age group, race, and PIR group as covariates. Matching uses a
# nearest-neighbor algorithm with a 1:1 ratio without replacement and a
# caliper of 0.2 standard deviations of the logit of the propensity score.
# Covariate balance before/after matching is assessed by standardized mean
# differences (SMD); an SMD <= 0.1 is considered indicative of adequate
# balance (see love_plot_smd.R). To avoid double adjustment of covariates
# already balanced through matching, the post-match logistic models are
# minimally adjusted, including only those covariates showing residual
# post-match imbalance (SMD > 0.1) other than the exposure (WWI) itself.
# ============================================================================

suppressMessages({ library(Matching); library(survey) })

## ---- Input: derived analytic cohort (n = 6,782; 208 cases) ------------------
## breast (0/1), Age_group (< 60 / >= 60), Race, PIR_group, Age (continuous),
## Education, Hypertension, Diabetes, WWI (continuous).
# load("analytic_cohort_main.RData")   # <- derived analytic file (placeholder)

## ---- Propensity score model: age group + race + PIR group -------------------
ps_model  <- glm(breast ~ Age_group + Race + PIR_group,
                 family = binomial(link = "logit"), data = cohort)
ps        <- predict(ps_model, type = "response")
logit_ps  <- log(ps / (1 - ps))

## ---- 1:1 nearest-neighbor matching without replacement ----------------------
## caliper = 0.2 * SD of the logit of the propensity score
set.seed(1)
m <- Match(Tr = cohort$breast, X = logit_ps, M = 1,
           caliper = 0.2 * sd(logit_ps),
           replace = FALSE, ties = TRUE,
           distance.tolerance = 1e-04, version = "fast")

matched <- cohort[unlist(m[c("index.treated", "index.control")]), ]
cat(sprintf("Matched cohort: %d participants (%d cases, %d controls)\n",
            nrow(matched),
            sum(matched$breast == 1), sum(matched$breast == 0)))

## ---- Post-match logistic models (minimal adjustment) -------------------------
## In this cohort the covariates with residual post-match imbalance
## (SMD > 0.1) other than WWI itself are: age (continuous), education,
## hypertension, and diabetes. These four are the only covariates entered
## beyond WWI in the adjusted post-match model.
matched$WWI_Q <- cut(matched$WWI, breaks = quantile(matched$WWI),
                     include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4"))

post_adj <- "Age + Education + Hypertension + Diabetes"   # minimal adjustment set
m1_cont <- glm(breast ~ WWI,     family = binomial(), data = matched)
m2_cont <- glm(as.formula(paste("breast ~ WWI +",  post_adj)), family = binomial(), data = matched)
m1_q    <- glm(breast ~ WWI_Q,   family = binomial(), data = matched)
m2_q    <- glm(as.formula(paste("breast ~ WWI_Q +", post_adj)), family = binomial(), data = matched)

or_ci <- function(fit, term) {
  est <- coef(fit)[term]
  ci  <- confint(fit, term)
  data.frame(OR = exp(est), lower = exp(ci[1]), upper = exp(ci[2]),
             P = summary(fit)$coef[term, "Pr(>|z|)"])
}
or_ci(m1_cont, "WWI");  or_ci(m2_cont, "WWI")
or_ci(m1_q, "WWI_QQ4"); or_ci(m2_q, "WWI_QQ4")

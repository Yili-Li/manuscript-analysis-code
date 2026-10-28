# ============================================================================
# Restricted cubic splines: WWI and breast cancer prevalence
# (manuscript Section 2.7)
#
# WWI modeled with 3 knots placed at the 10th, 50th, and 90th percentiles of
# its distribution (10.2, 11.3, 12.4 in the analytic cohort), in a
# multivariable logistic model adjusted for age, race, PIR, hypertension,
# diabetes, and CHD. The reference value (OR = 1) is set at the median WWI of
# the analytic cohort (11.3). Overall association and non-linearity tests are
# reported from the fitted spline. An exploratory inflection point is derived
# from the fitted curve; its 95% CI is estimated by nonparametric bootstrap
# with the percentile method (seed 20241115). Because the outer knot
# (90th percentile, WWI = 12.4) coincides numerically with the exploratory
# inflection estimate, the inflection is interpreted as a descriptive feature
# of the fitted curve rather than as a threshold.
# ============================================================================

suppressMessages({ library(rms); library(ggrcs); library(ggplot2) })

## ---- Input: derived analytic cohort (n = 6,782; 208 cases) ------------------
## breast (0/1), WWI (continuous), Age (years, continuous), Race, PIR,
## Hypertension, Diabetes, CHD.
# load("analytic_cohort_main.RData")   # <- derived analytic file (placeholder)

dd <- datadist(cohort)            # reference values default to medians
options(datadist = "dd")

fit <- lrm(breast ~ rcs(WWI, 3) + Age + Race + PIR +
             Hypertension + Diabetes + CHD,
           x = TRUE, y = TRUE, data = cohort)

anova(fit)        # WWI: overall association; nonlinearity test

## ---- OR curve across the observed WWI range, reference = median -------------
pdf("RCS_WWI_breast_cancer.pdf", width = 6, height = 4)
print(
  ggrcs(data = cohort, fit = fit, x = "WWI",
        xlab = "WWI (cm/sqrt(kg))", ylab = "OR (95% CI)") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme_bw()
)
dev.off()

## ---- Exploratory inflection: two-segment search + bootstrap percentile CI ---
covs <- "Age + Race + PIR + Hypertension + Diabetes + CHD"
two_segment_ll <- function(cut, d) {
  lo <- ifelse(d$WWI <= cut, d$WWI, cut)        # continuous two-piece linear
  hi <- ifelse(d$WWI >  cut, d$WWI - cut, 0)
  m  <- glm(as.formula(paste("breast ~ lo + hi +", covs)),
            data = d, family = binomial())
  as.numeric(logLik(m))
}
grid <- quantile(cohort$WWI, probs = seq(0.05, 0.95, by = 0.01))
point_estimate <- grid[which.max(sapply(grid, two_segment_ll, d = cohort))]

set.seed(20241115)
B <- 50
boot <- replicate(B, {
  d <- cohort[sample(nrow(cohort), replace = TRUE), ]
  if (length(unique(d$breast)) < 2) return(NA_real_)
  g <- quantile(d$WWI, probs = seq(0.05, 0.95, by = 0.01))
  g[which.max(sapply(g, two_segment_ll, d = d))]
})
ci <- quantile(boot, probs = c(0.025, 0.975), na.rm = TRUE)   # percentile CI
c(point = point_estimate, lower = ci[1], upper = ci[2])

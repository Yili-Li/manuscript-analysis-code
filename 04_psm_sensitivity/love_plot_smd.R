# ============================================================================
# Covariate balance before and after propensity-score matching (Love plot)
# (manuscript Section 2.5; Figure S1; Table S5)
#
# Variable-level standardized mean differences (SMD) are computed with
# tableone (categorical SMD follows Yang and Dalton 2012); an SMD <= 0.1 is
# considered indicative of adequate balance (Austin 2009). Educational
# attainment is categorized in 4 levels (less than 11th grade, combining
# "less than 9th grade" and "9th-11th grade"; high school graduate or GED;
# some college or associate degree; college graduate or above).
# ============================================================================

suppressMessages({ library(tableone); library(ggplot2) })

## ---- Input: pre-match and post-match cohort data frames ---------------------
## Both must contain the same baseline variables and a breast (0/1) column.
# load("analytic_cohort_main.RData");    pre  <- cohort    # n = 6,782
# load("analytic_cohort_matched.RData"); post <- matched   # n = 416 (1:1)

factor_vars <- c("Age_group", "Race", "Education", "PIR_group", "Marital",
                 "Alcohol", "Hypertension", "Diabetes", "CHD", "HeartFailure")
all_vars    <- c("Age_years", "WWI", factor_vars)

smd_at <- function(d) {
  tab <- CreateTableOne(vars = all_vars, factorVars = factor_vars,
                        strata = "breast", data = d,
                        test = FALSE, smd = TRUE)
  out <- ExtractSmd(tab)
  setNames(as.numeric(out[, 1]), rownames(out))
}

smd <- data.frame(Variable = all_vars,
                  SMD_pre  = as.numeric(smd_at(pre)[all_vars]),
                  SMD_post = as.numeric(smd_at(post)[all_vars]))
print(smd, digits = 3)

## ---- Love plot ----------------------------------------------------------------
plot_df <- data.frame(Variable = rep(smd$Variable, 2),
                      SMD   = c(smd$SMD_pre, smd$SMD_post),
                      Stage = rep(c("Before matching", "After matching"),
                                  each = nrow(smd)))

ggplot(plot_df, aes(x = SMD, y = reorder(Variable, abs(SMD)))) +
  geom_point(aes(color = Stage), size = 2) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey40") +
  labs(x = "Standardized mean difference", y = NULL, color = NULL) +
  theme_bw() + theme(panel.grid.minor = element_blank())
ggsave("FigureS1_love_plot.pdf", width = 7, height = 5)

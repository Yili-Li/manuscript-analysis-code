# ============================================================================
# Prespecified subgroup analyses, forest plot, and interaction tests
# (manuscript Sections 2.6 and 2.8; Figure 4; Table S4)
#
# Ten prespecified subgroup families: age (< 60 vs >= 60 years), race,
# educational attainment, PIR, marital status, alcohol use, hypertension,
# diabetes, CHD, and heart failure. Within-stratum odds ratios (WWI Q4 vs Q1)
# are crude (unadjusted) estimates. Effect modification is tested by adding a
# WWI x subgroup product term to an unadjusted logistic model; the family of
# ten interaction P-values is corrected for multiple testing using the
# Benjamini-Hochberg false-discovery rate (q < 0.05 considered suggestive).
# ============================================================================

suppressMessages({ library(jstable); library(ggplot2) })

## ---- Input: derived analytic cohort (n = 6,782; 208 cases) ------------------
## breast (0/1), WWI_Q (quartiles, Q1 reference), subgroup variables as factors.
# load("analytic_cohort_main.RData")   # <- derived analytic file (placeholder)

subgroup_vars <- c("Age_group", "Race", "Education", "PIR", "Marital",
                   "Alcohol", "Hypertension", "Diabetes", "CHD", "HeartFailure")

## ---- Crude within-stratum ORs (Q4 vs Q1) + interaction P-values -------------
res <- TableSubgroupMultiGLM(
  formula       = breast ~ WWI_Q,
  var_subgroups = subgroup_vars,
  var_cov       = NULL,
  family        = "binomial",
  data          = cohort)

print(res)

## ---- Benjamini-Hochberg FDR across the ten subgroup families ----------------
fam  <- res[is.na(res$Levels) | res$Levels == "", ]
q_bh <- p.adjust(as.numeric(fam$P.for.interaction), method = "BH")
data.frame(Subgroup = fam$Subgroup,
           P_interaction = fam$P.for.interaction,
           q_BH = q_bh)

## ---- Forest plot of within-stratum ORs --------------------------------------
lv <- res[!(is.na(res$Levels) | res$Levels == "") &
            !is.na(res$OR) & res$OR != "Reference", ]
lv$OR    <- as.numeric(lv$OR)
lv$Lower <- as.numeric(lv$Lower)
lv$Upper <- as.numeric(lv$Upper)
plot_df  <- data.frame(label = paste0(lv$Subgroup, ": ", lv$Levels),
                       OR = lv$OR, lower = lv$Lower, upper = lv$Upper)
plot_df$label <- factor(plot_df$label, levels = rev(plot_df$label))

ggplot(plot_df, aes(x = OR, y = label)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_x_log10() +
  labs(x = "OR (Q4 vs Q1, 95% CI)", y = NULL,
       title = "Subgroup analyses (crude)") +
  theme_bw() + theme(panel.grid.minor = element_blank())
ggsave("Figure4_subgroup_forest.pdf", width = 7, height = 6)

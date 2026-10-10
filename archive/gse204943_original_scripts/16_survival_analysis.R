# ==============================================================================
# Script 16: PMD Erosion Survival Analysis
# Project: OPMD Methylation Analysis (GSE204943)
# Author: Michael E. Troka
#
# Purpose:
#   Comprehensive survival analysis linking PMD erosion to clinical outcomes.
#   Creates pseudotime ordering, composite PMD erosion score (z-scored load +
#   inside-outside delta), Kaplan-Meier curves with multiple cutpoint
#   strategies (median/tertile/maxstat), univariate and multivariate Cox
#   regression adjusting for clinical variables, forest plots, and landmark
#   survival tables
#
# Inputs:
#   - pmd_load_df2
#   - pmd_sample_stats2
#   - survival metadata (Survival_in_months, Recurrence_status)
#   - clinical variables (Age, Sex, Stage, Nodal_status)
#
# Outputs:
#   - PMD_erosion_score_per_sample.csv
#   - pseudotime figures
#   - KM survival plots
#   - forest plots
#   - Table1_Univariate_Cox.docx
#   - Table2_Multivariate_Cox.docx
#   - Table3_KM_Landmark_Timepoints.docx
#
# Dependencies:
#   survival, survminer, dplyr, tidyr, ggplot2, scales, stringr, flextable
#
# Upstream: Script 10
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(stringr)
  library(survival)
  library(survminer)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

# Ensure PMD load df is clean
pmd_load_df2 <- pmd_load_df2 %>%
  mutate(
    Sample = as.character(Sample),
    Sample_Group = as.character(Sample_Group),
    Sample_Group = case_when(
      str_detect(Sample_Group, regex("tumor", ignore_case=TRUE))  ~ "1_Tumor",
      str_detect(Sample_Group, regex("opl", ignore_case=TRUE))    ~ "2_OPL",
      str_detect(Sample_Group, regex("normal|control", ignore_case=TRUE)) ~ "3_Normal",
      TRUE ~ Sample_Group
    ),
    Sample_Group = factor(Sample_Group, levels=c("3_Normal","2_OPL","1_Tumor"))
  )

# Join delta if useful later
pmd_all <- pmd_load_df2 %>%
  left_join(
    pmd_sample_stats2 %>% select(Sample, delta_in_minus_out),
    by="Sample"
  )




















# ----- pseudotime score -----
pmd_all <- pmd_all %>%
  mutate(
    PMD_pseudotime = scales::rescale(PMD_load_Mb, to=c(0,1))  # 0=normal-like, 1=most eroded
  ) %>%
  arrange(PMD_pseudotime)

# ----- plot trajectory with group colors -----
pal_groups <- c(
  "3_Normal"="#3C5488FF",
  "2_OPL"   ="#00A087FF",
  "1_Tumor" ="#E64B35FF"
)

p_pseudotime <- ggplot(pmd_all,
                       aes(x=PMD_pseudotime, y=PMD_load_Mb, color=Sample_Group)) +
  geom_point(size=2.2, alpha=0.9) +
  geom_smooth(method="loess", se=FALSE, linewidth=1) +
  scale_color_manual(values=pal_groups) +
  labs(
    title="Epigenomic pseudotime based on hypo-PMD erosion",
    x="PMD pseudotime (0 → 1)",
    y="Hypo-PMD coverage (Mb)",
    color=NULL
  ) +
  theme_classic(base_size=14)

ggsave(file.path(fig_dir,"Fig_pseudotime_PMD_load.pdf"),
       p_pseudotime, width=6.3, height=4.5, useDingbats=FALSE)

print(p_pseudotime)

# Optional: show how delta changes over pseudotime
p_pseudotime_delta <- ggplot(pmd_all,
                             aes(x=PMD_pseudotime, y=delta_in_minus_out, color=Sample_Group)) +
  geom_point(size=2.2, alpha=0.9) +
  geom_smooth(method="loess", se=FALSE, linewidth=1) +
  scale_color_manual(values=pal_groups) +
  labs(
    title="Inside–outside PMD methylation contrast across pseudotime",
    x="PMD pseudotime (0 → 1)",
    y=expression(Delta*beta~"(inside − outside)"),
    color=NULL
  ) +
  theme_classic(base_size=14)

ggsave(file.path(fig_dir,"Fig_pseudotime_delta.pdf"),
       p_pseudotime_delta, width=6.3, height=4.5, useDingbats=FALSE)

print(p_pseudotime_delta)




















# z-score both components and combine
pmd_all <- pmd_all %>%
  mutate(
    z_load  = as.numeric(scale(PMD_load_Mb)),
    z_delta = as.numeric(scale(-delta_in_minus_out)),  # negate so "more erosion" = higher
    PMD_erosion_score = (z_load + z_delta)/2
  )

# Plot score by group
p_score <- ggplot(pmd_all, aes(x=Sample_Group, y=PMD_erosion_score, fill=Sample_Group)) +
  geom_violin(trim=FALSE, alpha=0.55, color=NA) +
  geom_boxplot(width=0.18, outlier.shape=NA, alpha=0.95, linewidth=0.35) +
  geom_jitter(width=0.12, size=1.6, alpha=0.7) +
  scale_fill_manual(values=pal_groups) +
  labs(
    title="Composite PMD erosion score across biological groups",
    x=NULL, y="PMD erosion score (z-scaled)"
  ) +
  theme_classic(base_size=14) +
  theme(legend.position="none", axis.text.x=element_text(face="bold"))

ggsave(file.path(fig_dir,"Fig_PMD_erosion_score_by_group.pdf"),
       p_score, width=5.3, height=4.8, useDingbats=FALSE)

print(p_score)

# Save score table for downstream models
write.csv(pmd_all, file.path(results_dir,"PMD_erosion_score_per_sample.csv"),
          row.names=FALSE)















library(dplyr)
library(survival)
library(survminer)
library(stringr)

# Clean clinical variable whitespace before any filtering
pmd_all <- pmd_all %>%
  mutate(across(where(is.character), str_trim))

# CORRECTED event definition - the logic was inverted!
surv_df_tumor <- pmd_all %>%
  filter(Sample_Group == "1_Tumor", !is.na(Survival_in_months)) %>%
  mutate(
    # Fixed: "Recurrence" should be event=1, "No recurrence" should be event=0
    event = case_when(
      str_detect(Recurrence_status, "Recurrence") ~ 1,  # This is the recurrence
      Recurrence_status == "No recurrence" ~ 0,  # This is censored
      Recurrence_status == "LFU" ~ 0,  # Lost to follow-up = censored
      TRUE ~ 0
    )
  )

# Verify the corrected coding
cat("Event distribution:\n")
table(surv_df_tumor$Recurrence_status, surv_df_tumor$event)
cat("\nTotal samples:", nrow(surv_df_tumor),
    "\nEvents:", sum(surv_df_tumor$event),
    "\nCensored:", sum(surv_df_tumor$event == 0), "\n")

# Now continue with the analysis using corrected events
# 1. Median split
median_cut <- median(surv_df_tumor$PMD_erosion_score)
surv_df_tumor$erosion_median <- ifelse(surv_df_tumor$PMD_erosion_score >= median_cut,
                                       "High", "Low")

# 2. Tertile split (fixing the reference issue)
tertiles <- quantile(surv_df_tumor$PMD_erosion_score, c(0.33, 0.67))
surv_df_tumor$erosion_tertile <- case_when(
  surv_df_tumor$PMD_erosion_score < tertiles[1] ~ "Low",
  surv_df_tumor$PMD_erosion_score > tertiles[2] ~ "High",
  TRUE ~ "Medium"
)

# 3. Optimal cutpoint
library(maxstat)
maxstat_cut <- maxstat.test(
  Surv(Survival_in_months, event) ~ PMD_erosion_score,
  data = surv_df_tumor,
  smethod = "LogRank",
  pmethod = "HL",
  minprop = 0.20,
  maxprop = 0.80
)
cat("\nOptimal cutpoint:", maxstat_cut$estimate, "\n")
print(maxstat_cut)

surv_df_tumor$erosion_optimal <- ifelse(
  surv_df_tumor$PMD_erosion_score >= maxstat_cut$estimate,
  "High", "Low"
)

# Run survival analyses with CORRECTED event coding
km_median <- survfit(Surv(Survival_in_months, event) ~ erosion_median, data = surv_df_tumor)
km_optimal <- survfit(Surv(Survival_in_months, event) ~ erosion_optimal, data = surv_df_tumor)

# Log-rank tests
cat("\n=== Median split log-rank test ===\n")
print(survdiff(Surv(Survival_in_months, event) ~ erosion_median, data = surv_df_tumor))

cat("\n=== Optimal split log-rank test ===\n")
print(survdiff(Surv(Survival_in_months, event) ~ erosion_optimal, data = surv_df_tumor))

# Plot with median split
p_km_median <- ggsurvplot(
  km_median,
  data = surv_df_tumor,
  risk.table = TRUE,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = TRUE,
  xlab = "Time (months)",
  ylab = "Recurrence-free probability",
  title = paste0("RFS by PMD Erosion Score (Median split = ", round(median_cut, 3), ")"),
  legend.labs = c("Low PMD erosion", "High PMD erosion"),
  palette = c("#4575b4", "#d73027"),
  ggtheme = theme_classic(base_size = 14)
)
print(p_km_median)


# Cox regression with continuous PMD score
cox_continuous <- coxph(
  Surv(Survival_in_months, event) ~ PMD_erosion_score,
  data = surv_df_tumor
)
cat("\n=== Cox regression (continuous PMD score) ===\n")
summary(cox_continuous)

# Multivariate Cox with clinical variables
cox_multi <- coxph(
  Surv(Survival_in_months, event) ~ PMD_erosion_score + Age + Sex +
    Early_Late + Nodal_status,
  data = surv_df_tumor
)
cat("\n=== Multivariate Cox regression ===\n")
summary(cox_multi)




















# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(survminer)
  library(stringr)
  library(ggplot2)
})

# Calculate composite PMD erosion score
pmd_all <- pmd_all %>%
  mutate(
    z_load = as.numeric(scale(PMD_load_Mb)),
    z_delta = as.numeric(scale(-delta_in_minus_out)),
    PMD_erosion_score = (z_load + z_delta)/2
  )



# Create tumor-only survival dataset with correct event coding
surv_df_tumor <- pmd_all %>%
  filter(Sample_Group == "1_Tumor", !is.na(Survival_in_months)) %>%
  mutate(
    event = case_when(
      str_detect(Recurrence_status, "Recurrence") ~ 1,
      Recurrence_status == "No recurrence" ~ 0,
      Recurrence_status == "LFU" ~ 0,
      TRUE ~ 0
    )
  )

# Report sample characteristics
cat("Tumor samples:", nrow(surv_df_tumor),
    "\nEvents:", sum(surv_df_tumor$event),
    "\nCensored:", sum(surv_df_tumor$event == 0), "\n")



# Perform median split
median_cut <- median(surv_df_tumor$PMD_erosion_score)
surv_df_tumor$erosion_group <- factor(
  ifelse(surv_df_tumor$PMD_erosion_score >= median_cut,
         "High PMD erosion", "Low PMD erosion"),
  levels = c("Low PMD erosion", "High PMD erosion")
)

# Fit Kaplan-Meier model
km_fit <- survfit(Surv(Survival_in_months, event) ~ erosion_group,
                  data = surv_df_tumor)

# Log-rank test
survdiff(Surv(Survival_in_months, event) ~ erosion_group,
         data = surv_df_tumor)

# Generate publication figure
p_km <- ggsurvplot(
  km_fit,
  data = surv_df_tumor,
  risk.table = TRUE,
  pval = TRUE,
  pval.method = TRUE,
  conf.int = TRUE,
  xlab = "Time (months)",
  ylab = "Recurrence-free probability",
  title = "Recurrence-Free Survival by PMD Erosion Score in OSCC",
  subtitle = paste("Median split (cutoff =", round(median_cut, 3), ")"),
  legend.labs = c("Low PMD erosion", "High PMD erosion"),
  palette = c("#4575b4", "#d73027"),
  break.time.by = 20,
  risk.table.height = 0.25,
  ggtheme = theme_classic(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "top")
)

# Save figure
ggsave("Figure_RFS_PMD_erosion.pdf",
       plot = p_km$plot, width = 8, height = 6, dpi = 300)



# Univariate Cox regression (continuous PMD score)
cox_univariate <- coxph(
  Surv(Survival_in_months, event) ~ PMD_erosion_score,
  data = surv_df_tumor
)
summary(cox_univariate)

# Multivariate Cox regression
cox_multivariate <- coxph(
  Surv(Survival_in_months, event) ~ PMD_erosion_score + Age + Sex +
    Early_Late + Nodal_status,
  data = surv_df_tumor
)
summary(cox_multivariate)


# Extract hazard ratios and confidence intervals
cox_summary <- summary(cox_multivariate)

# PMD erosion score HR
pmd_hr <- exp(coef(cox_multivariate)["PMD_erosion_score"])
pmd_ci <- exp(confint(cox_multivariate)["PMD_erosion_score",])
pmd_p <- cox_summary$coefficients["PMD_erosion_score", "Pr(>|z|)"]

cat("PMD Erosion Score Results:\n",
    "HR =", round(pmd_hr, 2),
    "(95% CI:", round(pmd_ci[1], 2), "-", round(pmd_ci[2], 2), ")\n",
    "P-value =", format.pval(pmd_p, digits = 3), "\n")

# Model concordance
cat("Model C-index:", round(cox_summary$concordance[1], 3), "\n")




library(forcats)
library(dplyr)
library(survival)
library(survminer)

# 1. Clean data and RENAME the messy factor levels
surv_df_tumor_clean <- surv_df_tumor %>%
  filter(!is.na(Age), !is.na(Sex), !is.na(Early_Late), !is.na(Nodal_status)) %>%
  mutate(
    Sex = factor(Sex),
    Early_Late = factor(Early_Late),
    # RENAME LEVELS HERE: Remove spaces and special chars (+, _)
    Nodal_status = fct_recode(factor(Nodal_status),
                              "Node_Negative" = "Node _ve ",
                              "Node_Positive" = "Node +ve")
  ) %>%
  mutate_if(is.factor, fct_drop)

# 2. Refit model with the sanitized levels
cox_multi_clean <- coxph(
  Surv(Survival_in_months, event) ~ PMD_erosion_score + Age + Sex +
    Early_Late + Nodal_status,
  data = surv_df_tumor_clean
)

# 3. Plot - This should now work
ggforest(cox_multi_clean, data = surv_df_tumor_clean)


# install.packages("forestmodel")  # Moved to setup — install once before running
library(forestmodel)

forest_model(cox_multi_clean)




















# Load required packages
library(gtsummary)
library(flextable)
library(officer)

# Table 1: Univariate Cox Regression
tbl_univariate <- cox_univariate %>%
  tbl_regression(
    exponentiate = TRUE,
    label = list(PMD_erosion_score ~ "PMD Erosion Score"),
    estimate_fun = ~style_ratio(., digits = 2),
    pvalue_fun = ~style_pvalue(., digits = 3)
  ) %>%
  add_glance_source_note(
    include = c(n, nevent, concordance)
  ) %>%
  modify_header(label = "**Variable**") %>%
  modify_caption("**Table 1. Univariate Cox Regression Analysis**")

# Table 2: Multivariate Cox Regression
tbl_multivariate <- cox_multi_clean %>%
  tbl_regression(
    exponentiate = TRUE,
    label = list(
      PMD_erosion_score ~ "PMD Erosion Score",
      Age ~ "Age (years)",
      Sex ~ "Sex",
      Early_Late ~ "Stage",
      Nodal_status ~ "Nodal Status"
    ),
    estimate_fun = ~style_ratio(., digits = 2),
    pvalue_fun = ~style_pvalue(., digits = 3)
  ) %>%
  add_glance_source_note(
    include = c(n, nevent, concordance)
  ) %>%
  modify_header(label = "**Variable**") %>%
  modify_caption("**Table 2. Multivariate Cox Regression Analysis**")

# Save as PDF
library(webshot2)

# Convert to flextable for better PDF output
ft_univariate <- as_flex_table(tbl_univariate)
ft_multivariate <- as_flex_table(tbl_multivariate)

# Save tables to Word first, then convert to PDF
save_as_docx(ft_univariate, path = "Table1_Univariate_Cox.docx")
save_as_docx(ft_multivariate, path = "Table2_Multivariate_Cox.docx")

# Alternative: Create custom formatted table
library(kableExtra)

# Extract model results
extract_cox_results <- function(model, model_name) {
  summary_model <- summary(model)
  coef_df <- as.data.frame(summary_model$coefficients)
  conf_df <- as.data.frame(summary_model$conf.int)

  results <- data.frame(
    Variable = rownames(coef_df),
    HR = round(conf_df[,"exp(coef)"], 2),
    CI_lower = round(conf_df[,"lower .95"], 2),
    CI_upper = round(conf_df[,"upper .95"], 2),
    P_value = format.pval(coef_df[,"Pr(>|z|)"], digits = 3, eps = 0.001)
  )

  results$HR_CI <- paste0(results$HR, " (", results$CI_lower, "-", results$CI_upper, ")")
  results <- results[, c("Variable", "HR_CI", "P_value")]

  return(results)
}

# Create tables
table_uni <- extract_cox_results(cox_univariate, "Univariate")
table_multi <- extract_cox_results(cox_multi_clean, "Multivariate")

# Format for publication
library(knitr)
library(kableExtra)

# install.packages("tinytex")  # Moved to setup — install once before running
library(tinytex)

# Univariate table
kable(table_uni,
      col.names = c("Variable", "HR (95% CI)", "P-value"),
      caption = "Table 1. Univariate Cox Regression Analysis",
      format = "latex",
      booktabs = TRUE) %>%
  kable_styling(latex_options = c("hold_position")) %>%
  save_kable("Table1_Univariate.pdf")

# Multivariate table
kable(table_multi,
      col.names = c("Variable", "HR (95% CI)", "P-value"),
      caption = "Table 2. Multivariate Cox Regression Analysis",
      format = "latex",
      booktabs = TRUE) %>%
  kable_styling(latex_options = c("hold_position")) %>%
  footnote(general = paste("n =", cox_multi_clean$n,
                           "; Events =", cox_multi_clean$nevent,
                           "; C-index =", round(summary(cox_multi_clean)$concordance[1], 3))) %>%
  save_kable("Table2_Multivariate.pdf")



# Summary of KM survival object at clinically relevant time points
summary(km_median, times = c(12, 24, 36, 48, 60, 120))





















# -------------------------------------------------
# FINAL WORKING CODE – Publication-quality KM landmark table
# -------------------------------------------------
library(dplyr)
library(flextable)
library(officer)

# Time points
time_points <- c(12, 24, 36, 48, 60)

# Extract KM summary
km_sum <- summary(km_median, times = time_points)

# CORRECT INDEXING:
# km_sum stores results as: [High_12m, High_24m, High_36m, ..., Low_12m, Low_24m, Low_36m, ...]
# So: indices 1-5 are High group, indices 6-10 are Low group

n_times <- length(time_points)
high_idx <- 1:n_times
low_idx <- (n_times + 1):(n_times * 2)

# Build clean data frame
km_table <- data.frame(
  Time = paste0(time_points, " months"),

  High_n_risk = km_sum$n.risk[high_idx],
  High_events = km_sum$n.event[high_idx],
  High_RFS    = sprintf("%.1f%%", 100 * km_sum$surv[high_idx]),
  High_CI     = sprintf("(%.1f%%–%.1f%%)",
                        100 * km_sum$lower[high_idx],
                        100 * km_sum$upper[high_idx]),

  Low_n_risk = km_sum$n.risk[low_idx],
  Low_events = km_sum$n.event[low_idx],
  Low_RFS    = sprintf("%.1f%%", 100 * km_sum$surv[low_idx]),
  Low_CI     = sprintf("(%.1f%%–%.1f%%)",
                       100 * km_sum$lower[low_idx],
                       100 * km_sum$upper[low_idx])
)

# Verify: print to console to check against your summary output
print(km_table)

# Create flextable with CORRECT modern syntax
ft <- flextable(km_table) %>%
  set_header_labels(
    Time        = "Time",
    High_n_risk = "N at risk",
    High_events = "Events",
    High_RFS    = "RFS",
    High_CI     = "95% CI",
    Low_n_risk  = "N at risk",
    Low_events  = "Events",
    Low_RFS     = "RFS",
    Low_CI      = "95% CI"
  ) %>%
  add_header_row(
    values = c("", "High PMD Erosion (n=25)", "Low PMD Erosion (n=32)"),
    colwidths = c(1, 4, 4)
  ) %>%
  align(align = "center", part = "all") %>%
  bold(part = "header") %>%
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Arial", part = "all") %>%
  set_caption("Table 3. Recurrence-Free Survival at Landmark Time Points by PMD Erosion Score (Median Split)") %>%
  # CORRECT footnote syntax (new flextable)
  add_footer_lines(
    "RFS = recurrence-free survival; CI = confidence interval. Log-rank p = 0.022. Values are Kaplan–Meier estimates."
  ) %>%
  autofit()

# Preview
print(ft)

# Save as Word (perfect for submission)
save_as_docx(ft, path = "Table3_KM_Landmark_Timepoints.docx")

# Optional: Save as high-res PNG
save_as_image(ft, path = "Table3_KM_Landmark_Timepoints.png")

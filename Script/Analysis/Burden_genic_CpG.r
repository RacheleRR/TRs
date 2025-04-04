# Load libraries
library(broom)
library(dplyr)
library(gtsummary)
library(ggplot2)
library(readr)
library(data.table)
library(gprofiler2)
library(httr)
library(readxl)
library(stringr)

# Load data
ehdn_results_annotated <- read.delim("~/EHDN_DBSCAN_correct/Result/results_dbscan_with_QC_NOUHR/AFTER_DBSCAN/ehdn_DBSCAN_annotated.tsv")
ehdn_results_reordered <- read.delim("~/EHDN_DBSCAN_correct/Result/results_dbscan_with_QC_NOUHR/AFTER_DBSCAN/ehdn_DBSCAN_reorder.tsv")

# Merge and clean data
ehdn_results_annotated <- ehdn_results_annotated %>%
  left_join(ehdn_results_reordered %>% select(repeatID, ref, outlier_label, outlier_label2), by = "repeatID", suffix = c("", ".reordered")) %>%
  mutate(
    ref = coalesce(ref.reordered, ref),
    outlier_label = coalesce(outlier_label.reordered, outlier_label),
    outlier_label2 = coalesce(outlier_label2.reordered, outlier_label2)
  ) %>%
  select(-ref.reordered, -outlier_label.reordered, -outlier_label2.reordered) %>%
  mutate(count = sapply(strsplit(outliers, ";"), length))

ehdn_results_annotated <- ehdn_results_annotated %>%
  mutate(case_control = case_when(
    outlier_label == "Case" ~ 1,
    outlier_label == "Control" ~ 0,
    outlier_label == "mixed" ~ 2,
    TRUE ~ NA_real_
  ))


  ehdn_results_annotated <- ehdn_results_annotated %>% filter(count ==1 )



merge <- ehdn_results_annotated %>%
  filter(case_control == 1 | case_control == 0) %>%
  mutate(
    intergenic_genic = ifelse(region == "intergenic", "intergenic", "genic"),
    genic = ifelse(intergenic_genic == "genic", 1, 0),
    intergenic = ifelse(intergenic_genic == "intergenic", 1, 0),
    CpG = ifelse(grepl("CG|GC|CGC|CGG|CCG|GCC|GGC|GCG", motif), 1, 0),
    exonic = ifelse(region == "exonic", 1, 0),
    splicing = ifelse(region == "splicing", 1, 0),
    intronic = ifelse(region == "intronic", 1, 0),
    UTR = ifelse(region %in% c("UTR3", "UTR5"), 1, 0),
    UTR3 = ifelse(region == "UTR3", 1, 0),
    UTR5 = ifelse(region == "UTR5", 1, 0),
    upstream = ifelse(region %in% c("upstream", "upstream;downstream"), 1, 0),
    downstream = ifelse(region %in% c("downstream", "upstream;downstream"), 1, 0),
    ncRNA = ifelse(grepl("^ncRNA_", region), 1, 0)
  )

  merge <- merge %>%
  mutate(
    non_genic = ifelse(genic == 1, 0, 1),  # Inverse of genic
    non_CpG = ifelse(CpG == 1, 0, 1)       # Inverse of CpG
  )

  #burden done together
  model_full <- glm(
  case_control ~ genic + non_genic + CpG + non_CpG, 
  family = binomial, 
  data = merge
)
results <- tidy(model_full, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term %in% c("genic", "non_genic", "CpG", "non_CpG")) %>%
  mutate(
    term = factor(term, levels = c("genic", "non_genic", "CpG", "non_CpG")),
    is_sig = p.value < 0.05
  )

  ggplot(results, aes(x = term, y = estimate, fill = is_sig)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("gray", "#D6604D")) +
  labs(
    title = "Burden Analysis: Genic/Non-Genic and CpG/Non-CpG",
    x = "Category",
    y = "Odds Ratio (log scale)",
    fill = "Significant (p < 0.05)"
  ) +
  scale_y_log10() +  # Log scale for better visualization
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

  


# BURDEN analysis for genic
model_with_pred <- glm(case_control ~ genic, family = binomial(link = 'logit'), data = merge)
log_odds <- tidy(model_with_pred, conf.int = TRUE)
odds_ratio <- tidy(model_with_pred, exponentiate = TRUE, conf.int = TRUE)
tab_logistic <- bind_cols(log_odds, odds_ratio)
modified_tab_logistic <- tab_logistic %>%
  select(term...1, estimate...2, std.error...3, estimate...9, conf.low...13, conf.high...14, p.value...5) %>%
  rename(
    covariate = term...1,
    log_odds = estimate...2,
    SE = std.error...3,
    odds_ratio = estimate...9,
    lower_OR = conf.low...13,
    upper_OR = conf.high...14,
    p.val = p.value...5
  )
corrected_p_values <- p.adjust(modified_tab_logistic$p.val, method = "fdr")
modified_tab_logistic <- modified_tab_logistic %>%
  mutate(adj_p_val = corrected_p_values)

# Graph for genic
ggplot(modified_tab_logistic[-1, ], aes(x = covariate, y = odds_ratio, fill = adj_p_val < 0.05)) +
  geom_bar(stat = "identity", show.legend = FALSE, color = "black", width = 0.5) +
  geom_errorbar(aes(ymin = upper_OR, ymax = lower_OR), show.legend = FALSE, position = position_dodge(width = 0.5), size = 0.4, width = 0.2) +
  scale_fill_manual(values = c("#4393C3", "#D6604D")) +
  geom_hline(yintercept = 1, lty = 2) +
  theme_classic() +
  ylab("Odds Ratio") +
  xlab("") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    panel.border = element_rect(fill = NA)
  )

# BURDEN analysis for CpG
model_with_pred <- glm(case_control ~ CpG, family = binomial(link = 'logit'), data = merge)
log_odds <- tidy(model_with_pred, conf.int = TRUE)
odds_ratio <- tidy(model_with_pred, exponentiate = TRUE, conf.int = TRUE)
tab_logistic <- bind_cols(log_odds, odds_ratio)
modified_tab_logistic <- tab_logistic %>%
  select(term...1, estimate...2, std.error...3, estimate...9, conf.low...13, conf.high...14, p.value...5) %>%
  rename(
    covariate = term...1,
    log_odds = estimate...2,
    SE = std.error...3,
    odds_ratio = estimate...9,
    lower_OR = conf.low...13,
    upper_OR = conf.high...14,
    p.val = p.value...5
  )
corrected_p_values <- p.adjust(modified_tab_logistic$p.val, method = "fdr")
modified_tab_logistic <- modified_tab_logistic %>%
  mutate(adj_p_val = corrected_p_values)

# Graph for CpG
ggplot(modified_tab_logistic[-1, ], aes(x = covariate, y = odds_ratio, fill = adj_p_val < 0.05)) +
  geom_bar(stat = "identity", show.legend = FALSE, color = "black", width = 0.5) +
  geom_errorbar(aes(ymin = upper_OR, ymax = lower_OR), show.legend = FALSE, position = position_dodge(width = 0.5), size = 0.4, width = 0.2) +
  scale_fill_manual(values = c("#4393C3", "#D6604D")) +
  geom_hline(yintercept = 1, lty = 2) +
  theme_classic() +
  ylab("Odds Ratio") +
  xlab("") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    panel.border = element_rect(fill = NA)
)


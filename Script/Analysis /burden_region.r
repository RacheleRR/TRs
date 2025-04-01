#! Burden region 
#TODO :to understand how a TR in different parts of the genome may contribute to the risk of psychosis (quantify the association between the presence of a TR in different genomic regions and the likelihood of being a case)
#? clean and prepare dataframes
 #load libraries
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
 # load your data         
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
        outlier_label == "case" ~ 1,
        outlier_label == "control" ~ 0,
        outlier_label == "mixed" ~ 2,
        TRUE ~ NA_real_  # Handle any unexpected values
    ))

# Filter rows where case_control is either 1 or 0 and create a new data frame called merge
merge <- ehdn_results_annotated %>%
    filter(case_control == 1 | case_control == 0)

        merge <- merge%>% mutate(intergenic_genic = ifelse(region == "intergenic", "intergenic", "genic"))
        merge <- merge%>% mutate(genic = ifelse(intergenic_genic == "genic", 1,0))
        merge <- merge %>% mutate(intergenic = ifelse(intergenic_genic == "intergenic", 1,0))
        merge <- merge %>%  mutate(CpG = ifelse(grepl("CG|GC|CGC|CGG|CCG|GCC|GGC|GCG", motif), 1, 0))
        merge <- merge %>% mutate(exonic = ifelse(region == "exonic", 1,0))
        merge <- merge %>% mutate(splicing = ifelse(region == "splicing", 1,0))
        merge <- merge %>% mutate(intronic = ifelse(region == "intronic", 1,0))
        merge <- merge %>%  mutate(UTR = ifelse(region %in% c("UTR3", "UTR5"), 1, 0))
        merge <- merge %>%  mutate(UTR3 = ifelse(region =="UTR3", 1, 0))
        merge <- merge %>%  mutate(UTR5 = ifelse(region == "UTR5", 1, 0))
        merge <- merge %>%  mutate(upstream = ifelse(region %in% c("upstream", "upstream;downstream"), 1, 0))
        merge <- merge %>%  mutate(downstream = ifelse(region %in% c("downstream", "upstream;downstream"), 1, 0))
        merge <- merge %>%  mutate(ncRNA = ifelse(grepl("^ncRNA_", region), 1, 0))

#? LOGISTICAL REGRESSION ===  ALL TOGETHER (<- IS THE ONE I USED)

    model_with_pred <- glm(merge$case_control~merge$exonic + merge$intergenic +merge$UTR+ merge$upstream +merge$downstream+ merge$ncRNA+ merge$intronic,family= binomial(link = 'logit'),data = merge)

    # read data better
        log_odds <- tidy(model_with_pred, 
                        conf.int = TRUE)

        # convert log odds into odds ratio 
        odds_ratio <- tidy(model_with_pred,
                        exponentiate = TRUE,  
                        conf.int = TRUE)

        # unite them 
        tab_logistic <- bind_cols(log_odds, odds_ratio) 

        modified_tab_logistic <- tab_logistic %>%
        select(term...1, estimate...2, std.error...3, 
                estimate...9, conf.low...13, conf.high...14, p.value...5) %>%
        rename(covariate = term...1, 
                log_odds = estimate...2,
                SE = std.error...3,
                odds_ratio = estimate...9,
                lower_OR = conf.low...13, 
                upper_OR = conf.high...14,
                p.val = p.value...5)

        corrected_p_values <- p.adjust(modified_tab_logistic$p.val, method = "fdr")
        modified_tab_logistic <- modified_tab_logistic %>%
            mutate(adj_p_val = corrected_p_values)


#? GRAPH
   # Define the p-value for the "exonic" category

    ggplot(modified_tab_logistic[-1, ], aes(x = covariate, y = odds_ratio, fill = adj_p_val < 0.05)) +
      geom_bar(stat = "identity", show.legend = FALSE, color = "black", width = 0.5) +
      geom_errorbar(aes(ymin = upper_OR, ymax = lower_OR), show.legend = FALSE, 
                    position = position_dodge(width = 0.5), size = 0.4, width = 0.2) +
      scale_fill_manual(values = c("#4393C3", "#D6604D")) +
      geom_hline(yintercept = 1, lty = 2) +
      theme_classic() +
      ylab("Odds Ratio") +
      xlab("") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
            axis.text.y = element_text(size = 14), 
            axis.title.y = element_text(size = 14),
            panel.border = element_rect(fill = NA)) 

#            


#ONE BY ONE 
log_binominal <- function(Dep, Indep, Family, Data){
    model_with_pred <- glm(Dep ~ Indep, family = Family, data = Data)

    # read data better
    log_odds <- tidy(model_with_pred, conf.int = TRUE)

    # convert log odds into odds ratio 
    odds_ratio <- tidy(model_with_pred, exponentiate = TRUE, conf.int = TRUE)

    # unite them 
    tab_logistic <- bind_cols(log_odds, odds_ratio) 

    modified_tab_logistic <- tab_logistic %>%
        select(term...1, estimate...2, std.error...3, 
                estimate...9, conf.low...13, conf.high...14, p.value...5) %>%
        rename(covariate = term...1, 
                log_odds = estimate...2,
                SE = std.error...3,
                odds_ratio = estimate...9,
                lower_OR = conf.low...13, 
                upper_OR = conf.high...14,
                p.val = p.value...5)

    return(modified_tab_logistic)  
}
regions <- c("exonic", "intergenic", "UTR", "splicing", "intergenic", "genic")

# Use lapply to run the function for each region
results <- lapply(regions, function(region) {
  log_binominal(merge$Case_control, merge[[region]], Family = binomial(link = 'logit'), Data = merge)
})

# Optionally, combine results into a single data frame
combined_results <- bind_rows(results)



#EXTRA DONE LATER
# MULTIPLE REGIONS AT THE SAME TIME 
merge <- merge %>%
  mutate(region_type = case_when(
    region == "intergenic" ~ "intergenic",
    region == "intronic" ~ "intronic",
    region == "UTR" ~ "UTR",
    region == "exonic" ~ "exonic",
    region == "splicing" ~ "splicing",
    region == "ncRNA" ~ "ncRNA",
    region == "upstream" ~ "upstream",
    region == "downstream" ~ "downstream",
    TRUE ~ "other"
  ))

model <- glm(case_control ~ region_type, family = binomial(link = 'logit'), data = merge)
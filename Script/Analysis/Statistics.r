# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(kableExtra)
library(tidyverse)
library(readxl)
library(broom)
library(ggsignif)  # For adding significance annotations
library(patchwork)
library(FSA)
library(dunn.test)
library(gridExtra)  # For arranging multiple plots
library(purrr)
# LOAD dataframes
ehdn_results_annotated <- read.delim("~/EHDN_DBSCAN_correct/Result/results_dbscan_with_QC_NOUHR/AFTER_DBSCAN/ehdn_DBSCAN_annotated.tsv")
ehdn_results_reordered <- read.delim("~/EHDN_DBSCAN_correct/Result/results_dbscan_with_QC_NOUHR/AFTER_DBSCAN/ehdn_DBSCAN_reorder.tsv")
manifest_case_control_sex <- read.delim("/home/rachele/manifest_case_control_sex.csv")

# ADD LABELS FOR DETAILED PHENOTYPEM
    # Define input and manifest files
    input_file <- "ehdn_results_annotated"
    manifest_file <- "manifest_case_control_sex"

    # Load dataframes
    manifest_df <- get(manifest_file)
    input_df <- get(input_file)

    # Function to check the label for the outlier values
    check_outlier_label <- function(outlier_value, manifest_df) {
        outlier_values <- strsplit(outlier_value, ";")[[1]]
        labels <- sapply(outlier_values, function(val) {
            if (val %in% manifest_df$Sequencing_number) {
                return(manifest_df$Status[manifest_df$Sequencing_number == val])
            } else {
                return(NA)
            }
        })
        unique_labels <- unique(labels)
        if (length(unique_labels) > 1) {
            return("mixed")
        } else if (length(unique_labels) == 1) {
            return(unique_labels)
        } else {
            return(NA)
        }
    }

    # Function to get the labels (case or control) for the outlier values
    get_outlier_labels <- function(outlier_value, manifest_df) {
        outlier_values <- strsplit(outlier_value, ";")[[1]]
        labels <- sapply(outlier_values, function(val) {
            if (val %in% manifest_df$Sequencing_number) {
                return(manifest_df$Status[manifest_df$Sequencing_number == val])
            } else {
                return(NA)
            }
        })
        return(paste(labels, collapse = ", "))
    }

    # Apply functions to create new columns
    input_df$outlier_label_detail <- sapply(input_df$outliers, check_outlier_label, manifest_df = manifest_df)
    input_df$outlier_label_detail2 <- sapply(input_df$outliers, get_outlier_labels, manifest_df = manifest_df)
    assign(input_file, input_df)


    ehdn_results_annotated <- ehdn_results_annotated %>%
    left_join(ehdn_results_reordered %>% select(repeatID, ref, outlier_label, outlier_label2), by = "repeatID", suffix = c("", ".reordered")) %>%
    mutate(
        ref = coalesce(ref.reordered, ref),
        outlier_label = coalesce(outlier_label.reordered, outlier_label),
        outlier_label2 = coalesce(outlier_label2.reordered, outlier_label2)
    ) %>%
    select(-ref.reordered, -outlier_label.reordered, -outlier_label2.reordered) %>%
    mutate(count = sapply(strsplit(outliers, ";"), length))



# Preprocess
    # Add count columns
    ehdn_results_annotated <- ehdn_results_annotated %>%
        mutate(count_case = str_count(outlier_label2, "\\bCase\\b"),
            count_control = str_count(outlier_label2, "\\bControl\\b"),
            count_SCZ = str_count(outlier_label_detail2, "\\bFEP-SCZ\\b"),
            count_BD = str_count(outlier_label_detail2, "\\bFEP-BD\\b"),
            count_converter = str_count(outlier_label_detail2, "\\bConverter\\b"),
            count_non_converter = str_count(outlier_label_detail2, "\\bNon_Converter\\b"))

        ehdn_results_annotated <- ehdn_results_annotated%>%   mutate(CpG = ifelse(grepl("CG|GC|CGC|CGG|CCG|GCC|GGC|GCG", motif), 1, 0))
        

    # test
    # Summarize counts by region
    counts <- ehdn_results_annotated %>% group_by(region) %>%
        summarise(
            case_mix = sum(outlier_label %in% c("Case", "mixed")),
            control_mix = sum(outlier_label %in% c("Control", "mixed")),
            case_pure = sum(outlier_label == "Case"),
            control_pure = sum(outlier_label == "Control"),
            SCZ_mix = sum(outlier_label_detail %in% c("FEP-SCZ", "mixed")),
            BD_mix = sum(outlier_label_detail %in% c("FEP-BD", "mixed")),
            Converter_mix = sum(outlier_label_detail %in% c("Converter", "mixed")),
            Non_Converter_mix = sum(outlier_label_detail %in% c("Non_Converter", "mixed")),
            SCZ_pure = sum(outlier_label_detail == "FEP-SCZ"),
            BD_pure = sum(outlier_label_detail == "FEP-BD"),
            Converter_pure = sum(outlier_label_detail == "Converter"),
            Non_Converter_pure = sum(outlier_label_detail == "Non_Converter"),
            detail_case = sum(count_case),
            detail_control = sum(count_control),
            detail_SCZ = sum(count_SCZ),
            detail_BD = sum(count_BD),
            detail_converter = sum(count_converter),
            detail_non_converter = sum(count_non_converter)
        )


# FOR PSYCHO VS CONTROL 
    # Summarize case/control counts by region
    counts_case_control <- ehdn_results_annotated %>%
        group_by(region) %>%
        summarise(
            total = n(),
            case_mix = sum(outlier_label %in% c("Case", "mixed")),
            control_mix = sum(outlier_label %in% c("Control", "mixed")),
            case_pure = sum(outlier_label == "Case"),
            control_pure = sum(outlier_label == "Control"),
            detail_case = sum(count_case),
            detail_control = sum(count_control),
            CpG_control = sum(CpG[outlier_label == "Control"]),
            CpG_case = sum(CpG[outlier_label == "Case"]),
            Cpg = sum(CpG)
        )

    # Read manifest file and expand 'outliers' into separate rows
    manifest <- read.delim("~/EHDN_DBSCAN_correct/manifest_noUHRNA.tsv", header=FALSE) %>%
        rename(sample_id = V1, status = V2)

    unique_case_control_counts <- ehdn_results_annotated %>%
        separate_rows(outliers, sep = ";") %>%
        left_join(manifest, by = c("outliers" = "sample_id")) %>%
        group_by(region, status) %>%
        summarise(unique_samples = n_distinct(outliers), .groups = "drop") %>%
        pivot_wider(names_from = status, values_from = unique_samples, values_fill = list(unique_samples = 0)) %>%
        rename(total_unique_cases = Case, total_unique_controls = Control)

    # Merge with main summary counts
    counts_case_control <- counts_case_control %>%
        left_join(unique_case_control_counts, by = "region")

    # Calculate proportions and counts
    total_cases <- 302
    total_controls <- 75

    results <- counts_case_control %>%
        mutate(
            prop_pure_case = case_pure / total,
            prop_case_mix = case_mix / total,
            prop_pure_control = control_pure / total,
            prop_control_mix = control_mix / total,
            prop_cases_associated = detail_case / total_cases,
            prop_unique_cases_associated = total_unique_cases / total_cases,
            prop_controls_associated = detail_control / total_controls,
            prop_unique_controls_associated = total_unique_controls / total_controls,
            prop_case =detail_case/total ,
            prop_control = detail_control/total,
            prop_CpG_case = CpG_case/total,
            prop_CpG_control = CpG_control/total,
            prop_CpG = Cpg/total
        )






# FISHER TEST
# For TRs total 
    # Create a contingency table for each region for TR prop_BD_associated
    contingency_table_tr <- counts_case_control %>%
        select(
            region,
            case_pure,
            control_pure,
            total,  # Ensure total column is selected
            CpG_case,
            CpG_control
        ) %>%
        mutate(
            case_non_outlier = total - case_pure,  # Non-outlier TRs in cases
            control_non_outlier = total - control_pure,  # Non-outlier TRs in controls
            Non_CpG_case = case_pure - CpG_case,  # Non-CpG outlier TRs in cases
            Non_CpG_control = control_pure - CpG_control  # Non-CpG outlier TRs in controls
        )

    # View the contingency table
    print(contingency_table_tr)

# FOR CpG
    # Function to perform Fisher's exact test
    perform_fisher_test <- function(data) {
        fisher.test(matrix(c(data$CpG_case, data$CpG_control, data$Non_CpG_case, data$Non_CpG_control), nrow = 2))
    }

    # Apply the test to each region
    fisher_results_cpg <- contingency_table_tr %>%
        group_by(region) %>%
        group_modify(~ broom::tidy(perform_fisher_test(.x)))
    # Ungroup the dataframe
    fisher_results_cpg <- fisher_results_cpg %>% ungroup()

    # Adjust p-values for multiple testing
    fisher_results_cpg <- fisher_results_cpg %>%
        mutate(adjusted_p_value_bonferroni = p.adjust(p.value, method = "bonferroni")) %>%
        mutate(adjusted_p_value_fdr = p.adjust(p.value, method = "fdr"))


# FOR pure 
    # Function to perform Fisher's exact test
    perform_fisher_test <- function(data) {
        fisher.test(matrix(c(data$case_pure, data$control_pure, data$case_non_outlier, data$control_non_outlier), nrow = 2))
    }

    # Apply the test to each region
    fisher_results_tr <- contingency_table_tr %>%
        group_by(region) %>%
        group_modify(~ broom::tidy(perform_fisher_test(.x)))

    # Ungroup the dataframe
    fisher_results_tr <- fisher_results_tr %>% ungroup()

    # Adjust p-values for multiple testing
    fisher_results_tr <- fisher_results_tr %>%
    mutate(
        adjusted_p_value_fdr = p.adjust(p.value, method = "fdr"),  # Benjamini-Hochberg
        adjusted_p_value_bonferroni = p.adjust(p.value, method = "bonferroni")  # Bonferroni
    )




# DO IT FOR THEM with individuals 

    # Create a contingency table for each region
    contingency_tables <- counts_case_control %>%
        select(region, case_pure, control_pure) %>%
        mutate(
            case_non_outlier = total_cases - case_pure,
            control_non_outlier = total_controls - control_pure
        ) %>%
        select(region, case_pure, control_pure, case_non_outlier, control_non_outlier)

    perform_fisher_test <- function(data) {
        fisher.test(matrix(c(data$case_pure, data$control_pure, data$case_non_outlier, data$control_non_outlier), nrow = 2))
    }

    # Apply the test to each region
    fisher_results <- contingency_tables %>%
        group_by(region) %>%
        group_modify(~ broom::tidy(perform_fisher_test(.x)))

    # Ungroup the dataframe
    fisher_results <- fisher_results %>% ungroup()

    # Adjust p-values for multiple testing
    fisher_results <- fisher_results %>% mutate(adjusted_p_value_bonferroni = p.adjust(p.value, method = "bonferroni")) %>%
        mutate(adjusted_p_value_fdr = p.adjust(p.value, method = "fdr"))


# Mann witney 
# NON NORMALIZED 

    data_tr_case <- ehdn_results_annotated %>%
    select(region, count_case, count_control) %>%
    pivot_longer(
        cols = c(count_case, count_control),
        names_to = "group",
        values_to = "tr_count"
    )

    data_tr_case <- data_tr_case %>%
    mutate(group = as.factor(group)) %>% mutate(tr_count = as.numeric(tr_count)) # Convert group to factor


    # Perform Wilcoxon rank-sum test
    wilcoxon_results <-  wilcox.test(tr_count ~ group, data = data_tr_case)
    wilcoxon_results_df <- tibble(
    W = wilcoxon_results$statistic,
    p_value = wilcoxon_results$p.value
    )

    # Group by region and perform Wilcoxon rank-sum tests
    wilcox_results_per_region <- data_tr_case %>%
    group_by(region) %>%
    group_modify(~ {
        wilcox_result <- wilcox.test(tr_count ~ group, data = .x)
        return(tibble(
        W = wilcox_result$statistic,
        p.value = wilcox_result$p.value
        ))
    })
    # Ungroup the dataframe
    wilcox_results_per_region <- wilcox_results_per_region %>% ungroup()

    # Adjust p-values for multiple testing
    wilcox_results_per_region <- wilcox_results_per_region  %>% mutate(adjusted_p_value_bonferroni = p.adjust(p.value, method = "bonferroni")) %>%
        mutate(adjusted_p_value_fdr = p.adjust(p.value, method = "fdr"))

    print(wilcox_results_per_region)

    

# Normalized Analysis
    total_individuals <- c(
        case = 302,
        control = 76
    )

    data_long_normalized_case <- data_tr_case %>%
        mutate(
            group = gsub("count_", "", group),
            tr_per_individual = tr_count / total_individuals[group]
        )
    data_long_normalized_case <- data_long_normalized_case %>%  mutate(tr_per_individual = as.numeric(tr_per_individual))


    # Perform Wilcoxon rank-sum test on normalized data
    wilcoxon_results_normalized <-  wilcox.test(tr_per_individual ~ group, data = data_long_normalized_case)
    wilcoxon_results_normalized_df <- tibble(
        W = wilcoxon_results_normalized$statistic,
        p_value = wilcoxon_results_normalized$p.value
    )

    # Perform Wilcoxon rank-sum test


    # Group by region and perform Wilcoxon rank-sum tests on normalized data
    wilcox_results_per_region_normalized <- data_long_normalized_case %>%
    group_by(region) %>%
    group_modify(~ {
        data_filtered <- .x %>% filter(group %in% c("case", "control"))
        wilcox_result <- wilcox.test(tr_per_individual ~ group, data = data_filtered)
        return(tibble(
        W = wilcox_result$statistic,
        p.value = wilcox_result$p.value
        ))
    })
    # Ungroup the dataframe
    wilcox_results_per_region_normalized <- wilcox_results_per_region_normalized %>% ungroup()

    # Adjust p-values for multiple testing
    wilcox_results_per_region_normalized <- wilcox_results_per_region_normalized %>% mutate(adjusted_p_value_bonferroni = p.adjust(p.value, method = "bonferroni")) %>%
        mutate(adjusted_p_value_fdr = p.adjust(p.value, method = "fdr"))
   

    print(wilcox_results_per_region_normalized)




#
## NOW FOR DETAILED SCZ
    # Summarize SCZ counts by region
    counts_SCZ <- ehdn_results_annotated %>% group_by(region) %>%
        summarise(
            total = n(),
            SCZ_mix = sum(outlier_label_detail %in% c("FEP-SCZ", "mixed")),
            BD_mix = sum(outlier_label_detail %in% c("FEP-BD", "mixed")),
            Converter_mix = sum(outlier_label_detail %in% c("Converter", "mixed")),
            Non_Converter_mix = sum(outlier_label_detail %in% c("Non_Converter", "mixed")),
            SCZ_pure = sum(outlier_label_detail == "FEP-SCZ"),
            BD_pure = sum(outlier_label_detail == "FEP-BD"),
            Converter_pure = sum(outlier_label_detail == "Converter"),
            Non_Converter_pure = sum(outlier_label_detail == "Non_Converter"),
            detail_SCZ = sum(count_SCZ),
            detail_BD = sum(count_BD),
            detail_converter = sum(count_converter),
            detail_non_converter = sum(count_non_converter),
            CpG_SCZ = sum(CpG[outlier_label_detail == "FEP-SCZ"]),
            CpG_BD = sum(CpG[outlier_label_detail == "FEP-BD"]),
            CpG_Converter = sum(CpG[outlier_label_detail == "Converter"]),
            CpG_Non_Converter = sum(CpG[outlier_label_detail == "Non_Converter"]),
            CpG = sum(CpG)
        )

    # Read manifest file and expand 'outliers' into separate rows
    manifest <- read_excel("~/Downloads/Liste_samples_WGS_final.xlsx")

    unique_case_control_counts <- ehdn_results_annotated %>%
        separate_rows(outliers, sep = ";") %>%
        left_join(manifest, by = c("outliers" = "Sequencing_number")) %>%
        group_by(region, Status) %>%
        summarise(unique_samples = n_distinct(outliers), .groups = "drop") %>%
        pivot_wider(names_from = Status, values_from = unique_samples, values_fill = list(unique_samples = 0)) %>%
        rename(total_unique_SCZ = "FEP-SCZ", total_unique_BD = "FEP-BD", total_unique_Converter = Converter, total_unique_Non_Converter = Non_Converter)

    # Merge with main summary counts
    counts_SCZ <- counts_SCZ %>%
        left_join(unique_case_control_counts, by = "region")

    # Calculate proportions and counts
    total_SCZ <-  222
    total_BD<- 33
    total_converter <- 47
    total_non_converter <- 75

    # Calculate proportions and counts for SCZ
    results_SCZ <- counts_SCZ %>%
        mutate(
            prop_SCZ_pure = SCZ_pure / total,
            prop_SCZ_mix = SCZ_mix / total,
            prop_BD_pure = BD_pure / total,
            prop_BD_mix = BD_mix / total,
            prop_Converter_pure = Converter_pure / total,
            prop_Converter_mix = Converter_mix / total,
            prop_Non_Converter_pure = Non_Converter_pure / total,
            prop_Non_Converter_mix = Non_Converter_mix / total,
            prop_SCZ_associated = detail_SCZ / total_SCZ,
            prop_unique_SCZ_associated = total_unique_SCZ / total_SCZ,
            prop_BD_associated = detail_BD / total_BD,
            prop_unique_BD_associated = total_unique_BD / total_BD,
            prop_Converter_associated = detail_converter / total_converter,
            prop_unique_Converter_associated = total_unique_Converter / total_converter,
            prop_Non_Converter_associated = detail_non_converter / total_non_converter,
            prop_unique_Non_Converter_associated = total_unique_Non_Converter / total_non_converter,
            prop_CpG_SCZ = CpG_SCZ / total,
            prop_CpG_BD = CpG_BD / total,
            prop_CpG_Converter = CpG_Converter / total,
            prop_CpG_Non_Converter = CpG_Non_Converter / total,
            prop_CpG = CpG / total
        )

# FISHER 
# FOR THE ONCE WITH TRs
        pairwise_data <- results_SCZ %>%
            select(
                total,
                region,
                SCZ_pure,
                BD_pure,
                Converter_pure,
                Non_Converter_pure,
                CpG_SCZ,
                CpG_BD,
                CpG_Converter,
                CpG_Non_Converter
            ) %>%
            mutate(
                Non_CpG_SCZ = SCZ_pure - CpG_SCZ,  # Non-CpG outlier TRs in SCZ
                Non_CpG_BD = BD_pure - CpG_BD,      # Non-CpG outlier TRs in BD
                Non_CpG_Converter = Converter_pure - CpG_Converter,  # Non-CpG outlier TRs in Converter
                Non_CpG_Non_Converter = Non_Converter_pure - CpG_Non_Converter  # Non-CpG outlier TRs in Non_Converter
            )

    perform_pairwise_fisher <- function(data, groups, type = "pure") {
        results <- list()
        
        for (i in 1:(length(groups) - 1)) {
            for (j in (i + 1):length(groups)) {
                group1 <- groups[i]
                group2 <- groups[j]
                
                # Create a contingency table for the two groups
                if (type == "pure") {
                    contingency_table_pure <- data %>%
                        mutate(
                            group1_pure = .data[[paste0(group1, "_pure")]],
                            group2_pure = .data[[paste0(group2, "_pure")]],
                            group1_non_outlier = total - .data[[paste0(group1, "_pure")]],
                            group2_non_outlier = total - .data[[paste0(group2, "_pure")]]
                        ) %>%
                        select(
                            region,
                            group1_pure,
                            group2_pure,
                            group1_non_outlier,
                            group2_non_outlier
                        )
                    
                    # Perform Fisher's exact test for each region
                    fisher_results <- contingency_table_pure %>%
                        group_by(region) %>%
                        group_modify(~ {
                            # Construct the 2x2 contingency table
                            contingency_matrix <- matrix(c(
                                .x$group1_pure, .x$group2_pure,  # Group 1 and Group 2 counts
                                .x$group1_non_outlier, .x$group2_non_outlier  # Non-outlier counts
                            ), nrow = 2)
                            
                            # Perform Fisher's exact test
                            broom::tidy(fisher.test(contingency_matrix))
                        })
                    
                } else if (type == "CpG") {
                    contingency_table_CpG <- data %>%
                        mutate(
                            group1_CpG = .data[[paste0("CpG_", group1)]],
                            group2_CpG = .data[[paste0("CpG_", group2)]],
                            group1_Non_CpG = .data[[paste0("Non_CpG_", group1)]],
                            group2_Non_CpG = .data[[paste0("Non_CpG_", group2)]]
                        ) %>%
                        select(
                            region,
                            group1_CpG,
                            group2_CpG,
                            group1_Non_CpG,
                            group2_Non_CpG
                        )
                    
                    # Perform Fisher's exact test for each region
                    fisher_results <- contingency_table_CpG %>%
                        group_by(region) %>%
                        group_modify(~ {
                            # Construct the 2x2 contingency table
                            contingency_matrix <- matrix(c(
                                .x$group1_CpG, .x$group2_CpG,  # Group 1 and Group 2 CpG counts
                                .x$group1_Non_CpG, .x$group2_Non_CpG  # Non-CpG counts
                            ), nrow = 2)
                            
                            # Perform Fisher's exact test
                            broom::tidy(fisher.test(contingency_matrix))
                        })
                }
                
                # Store the results
                results[[paste(group1, "vs", group2)]] <- fisher_results
            }
        }
        
        return(results)
    }

    # Define the groups to compare
    groups_to_compare <- c("SCZ", "BD", "Converter", "Non_Converter")

    # Perform pairwise Fisher's exact tests for pure TRs
    pairwise_results_pure <- perform_pairwise_fisher(pairwise_data, groups_to_compare, type = "pure")

    # Perform pairwise Fisher's exact tests for CpG-related TRs
    pairwise_results_CpG <- perform_pairwise_fisher(pairwise_data, groups_to_compare, type = "CpG")


    # Combine results for pure TRs into a single dataframe
    summary_table_pure <- do.call(rbind, lapply(names(pairwise_results_pure), function(comparison) {
        pairwise_results_pure[[comparison]] %>%
            mutate(Comparison = comparison)
    }))

    # Combine results for CpG-related TRs into a single dataframe
    summary_table_CpG <- do.call(rbind, lapply(names(pairwise_results_CpG), function(comparison) {
        pairwise_results_CpG[[comparison]] %>%
            mutate(Comparison = comparison)
    }))

    # Adjust p-values across all comparisons for pure TRs
    summary_table_pure <- summary_table_pure %>%
        mutate(
            adjusted_p_value_fdr = p.adjust(p.value, method = "fdr"),  # Benjamini-Hochberg (FDR)
            adjusted_p_value_bonferroni = p.adjust(p.value, method = "bonferroni")  # Bonferroni (FWER)
        )

    # Adjust p-values across all comparisons for CpG-related TRs
    summary_table_CpG <- summary_table_CpG %>%
        mutate(
            adjusted_p_value_fdr = p.adjust(p.value, method = "fdr"),  # Benjamini-Hochberg (FDR)
            adjusted_p_value_bonferroni = p.adjust(p.value, method = "bonferroni")  # Bonferroni (FWER)
        )

# 
# NOW DO IT WITH INDIVIDUALS 
    pairwise_data <- counts_SCZ %>%select(
            region,
            SCZ_pure,
            BD_pure,
            Converter_pure,
            Non_Converter_pure) %>%
            mutate(
                SCZ_non_outlier = total_SCZ - SCZ_pure,  # Non-outlier individuals in SCZ
                BD_non_outlier = total_BD - BD_pure,      # Non-outlier individuals in BD
                Converter_non_outlier = total_converter - Converter_pure,  # Non-outlier individuals in Converter
                Non_Converter_non_outlier = total_non_converter - Non_Converter_pure  # Non-outlier individuals in Non_Converter
            )
    


    # Function to perform pairwise Fisher's exact tests
        perform_pairwise_fisher_individuals <- function(data, groups) {
            results <- list()
            
            for (i in 1:(length(groups) - 1)) {
                for (j in (i + 1):length(groups)) {
                    group1 <- groups[i]
                    group2 <- groups[j]
                    
                    # Create a contingency table for the two groups
                    contingency_table <- data %>%
                        select(
                            region,
                            group1_pure = paste0(group1, "_pure"),
                            group2_pure = paste0(group2, "_pure"),
                            group1_non_outlier = paste0(group1, "_non_outlier"),
                            group2_non_outlier = paste0(group2, "_non_outlier")
                        )
                    
                    # Perform Fisher's exact test for each region
                    fisher_results <- contingency_table %>%
                        group_by(region) %>%
                        group_modify(~ broom::tidy(fisher.test(matrix(c(
                            .x$group1_pure, .x$group2_pure,  # Group 1 and Group 2 counts
                            .x$group1_non_outlier, .x$group2_non_outlier  # Non-outlier counts
                        ), nrow = 2))))
                        
                    # Store the results
                    results[[paste(group1, "vs", group2)]] <- fisher_results
                }
            }
            
            return(results)
        }


        # Define the groups to compare
        groups_to_compare <- c("SCZ", "BD", "Converter", "Non_Converter")

        # Perform pairwise Fisher's exact tests
        pairwise_results <- perform_pairwise_fisher_individuals(pairwise_data, groups_to_compare)

        summary_table_FISHER_SCZ <- do.call(rbind, lapply(names(pairwise_results), function(comparison) {
            pairwise_results[[comparison]] %>%
                mutate(Comparison = comparison)
        }))

        summary_table_FISHER_SCZ <- summary_table_FISHER_SCZ %>%
            mutate(
                adjusted_p_value_fdr = p.adjust(p.value, method = "fdr"),
                adjusted_p_value_bonferroni = p.adjust(p.value, method = "bonferroni")
            )



# 
# Kruskal 
# non normalized
    # Prepare data
    # Reshape the data
    data_tr <- ehdn_results_annotated %>%
        select(region, count_SCZ, count_BD,count_converter,count_non_converter) %>%
        pivot_longer(
            cols = c(count_SCZ, count_BD,count_converter,count_non_converter),
            names_to = "group",
            values_to = "tr_count"
        )
# manually 
    # # Perform Kruskal-Wallis test
    # kruskal_result <- kruskal.test(tr_count ~ group, data = data_tr)
    # print(kruskal_result)

    # # If significant, perform post-hoc Dunn's test
    # if (kruskal_result$p.value < 0.05) {
    #     library(dunn.test)
    #     dunn_result <- dunn.test(data_tr$tr_count, data_tr$group, method = "bh")
    #     print(dunn_result)
    # }

    # # Group by region and perform Kruskal-Wallis tests
    # krusk_results_per_region <- data_tr %>%
    #     group_by(region) %>%
    #     group_modify(~ {
    #         kruskal_result_region <- kruskal.test(tr_count ~ group, data = .x)
    #         return(tibble(statistic = kruskal_result_region$statistic, p.value = kruskal_result_region$p.value, method = "Kruskal-Wallis"))
    #     })

    # #if significant perform post-hoc Dunn's test
    # krusk_results_per_region <- data_tr %>%
    #     group_by(region) %>%
    #     group_modify(~ {
    #         kruskal_result_region <- kruskal.test(tr_count ~ group, data = .x)
    #         if (kruskal_result_region$p.value < 0.05) {
    #             dunn_result_region <- dunn.test(.x$tr_count, .x$group, method = "bh")
    #             return(tibble(statistic = dunn_result_region$statistic, p.value = dunn_result_region$p.value, method = "Dunn"))
    #         } else {
    #             return(tibble(statistic = kruskal_result_region$statistic, p.value = kruskal_result_region$p.value, method = "Kruskal-Wallis"))
    #         }
    #     })

    # # View results
    # print(krusk_results_per_region)

# NORMLIZE 
    data_tr <- ehdn_results_annotated %>%
        select(region, count_SCZ, count_BD,count_converter,count_non_converter) %>%
        pivot_longer(
            cols = c(count_SCZ, count_BD,count_converter,count_non_converter),
            names_to = "group",
            values_to = "tr_count"
        )


    # Define total individuals per group
    total_individuals <- c(
        SCZ = 229,
        BD = 33,
        converter = 50,
        nonConverter = 78
    )


    data_long <- data_tr %>%
        mutate(group = gsub("count_", "", group)) %>%
        mutate( group =gsub("non_converter", "nonConverter",group))
  

    # Normalize TR counts
    data_long_normalized <- data_long %>%
        mutate(tr_per_individual = tr_count / total_individuals[group])

# manually 
    # Perform Kruskal-Wallis test on normalized data
    # kruskal_result <- kruskal.test(tr_per_individual ~ group, data = data_long_normalized)
    # print(kruskal_result)
    # If significant, perform post-hoc Dunn's test
    # if (kruskal_result$p.value < 0.05) {
    #     library(dunn.test)
    #     dunn_result <- dunn.test(data_long_normalized$tr_per_individual, data_long_normalized$group, method = "bh")
    #     print(dunn_result)
    # }

    # # Group by region and perform Kruskal-Wallis tests on normalized data
    # krusk_results_per_region_normalized <- data_long_normalized %>%
    # group_by(region) %>%
    # group_modify(~ {
    #     kruskal_result <- kruskal.test(tr_per_individual ~ group, data = .x)
    #     return(tibble(statistic = kruskal_result$statistic, p.value = kruskal_result$p.value, method = "Kruskal-Wallis"))
    # })
    # # Ungroup the dataframe
    # krusk_results_per_region_normalized <- krusk_results_per_region_normalized %>% ungroup()

    # #adjust p-values for multiple testing
    # krusk_results_per_region_normalized <- krusk_results_per_region_normalized %>% mutate(adjusted_p_value_bonferroni = p.adjust(p.value, method = "bonferroni")) %>%
    #     mutate(adjusted_p_value_fdr = p.adjust(p.value, method = "fdr"))

    # # View results
    # print(krusk_results_per_region_normalized) 


# AUTOMATICESD
    # Function for Kruskal-Wallis and Dunn's tests
    perform_kruskal_dunn_test <- function(data, response_col, group_col) {
    kruskal_formula <- as.formula(paste(response_col, "~", group_col))
    kruskal_result <- kruskal.test(kruskal_formula, data = data)
    print(kruskal_result)
    
    dunn_test_result <- NULL
    if (kruskal_result$p.value <= 0.05) {
        dunn_test_result <- dunn.test(data[[response_col]], data[[group_col]], method = "bh",altp = TRUE)
        print(dunn_test_result)
    }
    
    list(kruskal_result = kruskal_result, dunn_test_result = dunn_test_result)
    }


    # Function for Kruskal-Wallis and Dunn's tests by region
    perform_kruskal_dunn_test_region <- function(data, response_col, group_col) {
    results <- list()
    
    for (region in unique(data$region)) {
        region_data <- data %>% filter(region == !!region)
        kruskal_formula <- as.formula(paste(response_col, "~", group_col))
        kruskal_result <- kruskal.test(kruskal_formula, data = region_data)
        print(kruskal_result)
        
        dunn_test_result <- NULL
        if (kruskal_result$p.value <= 0.05) {
        dunn_test_result <- dunn.test(region_data[[response_col]], region_data[[group_col]], method = "bh", altp = TRUE)
        print(dunn_test_result)
        }
        
        results[[region]] <- list(kruskal_result = kruskal_result, dunn_test_result = dunn_test_result)
    }
    
    return(results)
    }

# execution 
    # Non-normalized data
    kruskal_result <- perform_kruskal_dunn_test(data_tr, "tr_count", "group")
    krusk_results_per_region <- perform_kruskal_dunn_test_region(data_tr, "tr_count", "group")

    # Normalize data
    kruskal_result_normalized <- perform_kruskal_dunn_test(data_long_normalized, "tr_per_individual", "group")
    krusk_results_per_region_normalized <- perform_kruskal_dunn_test_region(data_long_normalized, "tr_per_individual", "group")



#
#
# VISULAIZATION
# PHYSO VS CONTROL
# FISHER
    # Combine TR and CpG results into one dataframe
    combined_results <- bind_rows(
    fisher_results_tr %>% mutate(test_type = "TR Outliers"),
    fisher_results_cpg %>% mutate(test_type = "CpG-Associated Outliers")
    )

    # Create Manhattan-style plot
    manhattan_plot <- ggplot(combined_results, aes(x = region, y = -log10(adjusted_p_value_fdr))) +
    geom_point(aes(color = adjusted_p_value_fdr < 0.05), size = 3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    facet_wrap(~test_type, ncol = 1) +  # Separate TR and CpG into panels
    scale_color_manual(values = c("gray", "blue"), name = "Significant (FDR < 0.05)") +
    labs(
        title = "Adjusted P-values for TR Outliers vs. CpG-Associated Outliers",
        x = "Region",
        y = "-log10(FDR-adjusted p-value)"
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold")
    )

    print(manhattan_plot)

    # Volcano plot for TR outliers
    volcano_tr <- ggplot(fisher_results_tr, aes(x = log2(estimate), y = -log10(adjusted_p_value_fdr))) +
    geom_point(aes(color = adjusted_p_value_fdr < 0.05), alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_color_manual(values = c("gray", "blue")) +
    labs(title = "TR Outliers", x = "log2(Odds Ratio)", y = "-log10(FDR p-value)") +
    theme_minimal()

    # Volcano plot for CpG-associated outliers
    volcano_cpg <- ggplot(fisher_results_cpg, aes(x = log2(estimate), y = -log10(adjusted_p_value_fdr))) +
    geom_point(aes(color = adjusted_p_value_fdr < 0.05), alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_color_manual(values = c("gray", "blue")) +
    labs(title = "CpG-Associated Outliers", x = "log2(Odds Ratio)", y = "") +
    theme_minimal()

    # Arrange side-by-side
    combined_volcano <- volcano_tr + volcano_cpg + plot_layout(guides = "collect")
    print(combined_volcano)

    volcano_tr + 
    ggrepel::geom_text_repel(
        data = filter(combined_results, adjusted_p_value_fdr < 0.05),
        aes(label = region), 
        box.padding = 0.5, 
        max.overlaps = 20
    )


    # Dot plot 
    dot_plot <- ggplot(combined_results, aes(x = -log10(adjusted_p_value_fdr), y = region)) +
    geom_point(aes(color = test_type, shape = adjusted_p_value_fdr < 0.05), size = 3) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_color_manual(values = c("CpG-Associated Outliers" = "darkgreen", "TR Outliers" = "purple")) +
    scale_shape_manual(values = c(16, 17)) +  # Shapes for non-sig vs. sig
    labs(
        title = "Significance by Region and Test Type",
        x = "-log10(FDR-adjusted p-value)",
        y = "Region",
        color = "Test Type",
        shape = "Significant (FDR < 0.05)"
    ) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8))

    print(dot_plot)



# FOR proportions
# for pure case and control
    # Filter and reshape data
    stacked_pure_data <- results %>%
    select(region, prop_pure_case, prop_pure_control) %>%
    pivot_longer(
        cols = starts_with("prop_"),
        names_to = "category",
        values_to = "proportion"
    )

    # Define a pretty color palette
    pure_colors <- c("#66C2A5", "#FC8D62")  # Teal and orange

    # Plot
    stacked_bar_pure <- ggplot(stacked_pure_data, aes(x = region, y = proportion, fill = category)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = pure_colors, labels = c("Pure Case", "Pure Control")) +  # Use `scale_fill_manual` for discrete data
    labs(
        title = "Proportions of Pure Cases and Controls by Region",
        x = "Region",
        y = "Proportion",
        fill = "Category"
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
    )

    print(stacked_bar_pure)

# for CpG_
    # Filter and reshape data
    stacked_cpg_data <- results %>%
    select(region, prop_CpG_case, prop_CpG_control) %>%
    pivot_longer(
        cols = starts_with("prop_"),
        names_to = "category",
        values_to = "proportion"
    )

    # Define a custom color palette
    cpg_colors <- c("#8DA0CB", "#E78AC3")  # Purple and pink

    # Plot
    stacked_bar_cpg <- ggplot(stacked_cpg_data, aes(x = region, y = proportion, fill = category)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = cpg_colors, labels = c("CpG Case", "CpG Control")) +  # Use `scale_fill_manual` for discrete data
    labs(
        title = "Proportions of CpG-Associated Cases and Controls by Region",
        x = "Region",
        y = "Proportion",
        fill = "Category"
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
    )

    print(stacked_bar_cpg)

# all proportions
    # Filter and reshape data
    heatmap_data <- results %>%
    select(region, prop_pure_case, prop_pure_control, prop_CpG_case, prop_CpG_control) %>%
    pivot_longer(
        cols = starts_with("prop_"),
        names_to = "category",
        values_to = "proportion"
    )

    # Plot
    heatmap <- ggplot(heatmap_data, aes(x = region, y = category, fill = proportion)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(
        title = "Proportions by Region and Category",
        x = "Region",
        y = "Category",
        fill = "Proportion"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

    print(heatmap)



#
# FOR SCZ 
# FISHER
# dot_plot
    dot_plot_SCZ <- ggplot(summary_table_pure, aes(x = -log10(adjusted_p_value_fdr), y = region)) +
        geom_point(aes(color = Comparison , shape = adjusted_p_value_fdr < 0.05), size = 3) +
        geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
        
        scale_shape_manual(values = c(16, 17)) +  # Shapes for non-sig vs. sig
        labs(
            title = "Significance by Region and Test Type",
            x = "-log10(FDR-adjusted p-value)",
            y = "Region",
            color = "Test Type",
            shape = "Significant (FDR < 0.05)"
        ) +
        theme_bw() +
        theme(axis.text.y = element_text(size = 8))

    print(dot_plot_SCZ)

# Volcano plot for TR outliers
    combined_results <- bind_rows(
    summary_table_pure %>% mutate(Test_Type = "Pure TRs"),
    summary_table_CpG %>% mutate(Test_Type = "CpG TRs"),
    summary_table_FISHER_SCZ %>% mutate(Test_Type = "Individuals")
    ) %>%
    mutate(
        log_OR = log2(estimate),
        log_p = -log10(adjusted_p_value_fdr),
        Significance = ifelse(adjusted_p_value_fdr < 0.05, "FDR < 0.05", "NS")
    )

    vulcano_plot_SCZ <- ggplot(combined_results, aes(x = log_OR, y = log_p, color = Significance)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.8, size = 3) +
    scale_color_manual(values = c("FDR < 0.05" = "red", "NS" = "grey60")) +
    facet_grid(Test_Type ~ Comparison, scales = "free") +
    theme_bw() +
    labs(
        x = "Log2(Odds Ratio)", 
        y = "-Log10(FDR-adjusted p)",
        title = "Volcano Plots of Pairwise Comparisons"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
vulcano_plot_SCZ +
    ggrepel::geom_text_repel(
        data = filter(combined_results, Significance == "FDR < 0.05"),  # Use the same significance threshold as the color aesthetic
        aes(label = region),
        box.padding = 0.5,   # Adjust spacing around labels
        max.overlaps = 20,   # Increase if you want more labels to show
        size = 3,            # Adjust font size
        color = "black"      # Ensure labels are visible against red points
    )
# Heatmap of log odds ratios + p value
    plot_heatmap <- function(test_type) {
    df <- combined_results %>% filter(Test_Type == test_type)
    
    ggplot(df, aes(x = Comparison, y = region, fill = log_OR)) +
        geom_tile(color = "white") +
        geom_text(aes(label = ifelse(Significance == "FDR < 0.05", "*", "")), 
                vjust = 0.8, size = 6) +
        scale_fill_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = 0, name = "Log2(OR)"
        ) +
        ggtitle(test_type) +
        theme_minimal() +
        theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")
        )
    }

    # Generate individual plots
    pure_plot <- plot_heatmap("Pure TRs")
    cpg_plot <- plot_heatmap("CpG TRs")
    individual_plot <- plot_heatmap("Individuals")

    # Arrange using patchwork
    library(patchwork)
    heatmap_combined <- (pure_plot / cpg_plot / individual_plot) + 
    plot_layout(guides = "collect")

# lollipop plot for effect size 
    lolli_SCZ <-ggplot(combined_results, aes(x = log_OR, y = reorder(region, log_OR))) +
    geom_segment(aes(xend = 0, yend = region), color = "grey50") +
    geom_point(aes(
        color = ifelse(adjusted_p_value_fdr < 0.05, "Sig", "Non-Sig"),
        size = 3
    )) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("Sig" = "red", "Non-Sig" = "grey60")) +
    facet_grid(Test_Type ~ Comparison, scales = "free_x") +
    labs(
        x = "Log2(Odds Ratio)", 
        y = "Genomic Region",
        color = "Significance"
    ) +
    theme_bw() +
    theme(
        axis.text.y = element_text(size = 8),
        strip.text = element_text(face = "bold")
    )


 # EXPLAINTION
# Positive Effect (log(OR) > 0):

  #      The first group in the comparison (e.g., SCZ) has higher odds of the outcome (e.g., having outlier TRs) compared to the second group (e.g., BD).

   #     Example: If SCZ vs. BD has a positive log(OR), it means SCZ individuals are more likely to have outlier TRs than BD individuals.

   # Negative Effect (log(OR) < 0):

 #       The first group in the comparison has lower odds of the outcome compared to the second group.

   #     Example: If SCZ vs. BD has a negative log(OR), it means SCZ individuals are less likely to have outlier TRs than BD individuals. --> -->

# dot_plot  for size effect and Direction
    plot_data <- combined_results %>%
    mutate(
        Comparison = factor(Comparison),
        Direction = ifelse(log_OR > 0, "Positive", "Negative")
    )

    dot_plot_size_effect_SCZ <-ggplot(plot_data, aes(x = Comparison, y = region)) +
    geom_point(aes(
        size = -log10(adjusted_p_value_fdr), 
        color = Direction,
        alpha = ifelse(Significance == "FDR < 0.05", 1, 0.3)
    )) +
    scale_size_continuous(range = c(2, 8)) +
    scale_color_manual(values = c("Positive" = "#d73027", "Negative" = "#4575b4")) +
    scale_alpha_identity() +
    facet_grid(Test_Type ~ ., scales = "free_y", space = "free") +
    theme_minimal() +
    labs(
        x = "Pairwise Comparison",
        y = "Genomic Region",
        size = "-Log10(p)",
        color = "Effect Direction"
    ) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(color = "grey90")
    )


# FOR PROPORTIONS
# pure SCZ
    stacked_bar_SCZ <- results_SCZ %>%
    select(region, prop_SCZ_pure, prop_BD_pure, prop_Converter_pure, prop_Non_Converter_pure) %>%
    pivot_longer(
    cols = starts_with("prop_"),
    names_to = "category",
    values_to = "proportion"
    )
    # Define a custom color palette
    scz_colors <- c("#8DA0CB", "#E78AC3", "#66B2FF", "#FFCC66")  # Purple, pink, blue, orange
    # Plot
    stacked_bar_SCZ <- ggplot(stacked_bar_SCZ, aes(x = region, y = proportion, fill = category)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = scz_colors, labels = c("SCZ", "BD", "Converter", "Non-Converter")) +  # Use `scale_fill_manual` for discrete data
    labs(
    title = "Proportions of SCZ-Associated Cases by Region",
    x = "Region",
    y = "Proportion",
    fill = "Category"
    ) +
    theme_bw() +
    theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
    )
    print(stacked_bar_SCZ)

# for cpg scz

    stacked_bar_SCZ_CpG <- results_SCZ %>%
    select(region, prop_CpG_SCZ, prop_CpG_BD, prop_CpG_Converter, prop_CpG_Non_Converter) %>%
    pivot_longer(
    cols = starts_with("prop_CpG_"),
    names_to = "category",
    values_to = "proportion"
    )   
    # Define a custom color palette
    scz_colors <- c("#8DA0CB", "#E78AC3", "#66B2FF", "#FFCC66")  # Purple, pink, blue, orange
    # Plot
    stacked_bar_SCZ_CpG <- ggplot(stacked_bar_SCZ_CpG, aes(x = region, y = proportion, fill = category)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = scz_colors, labels = c("SCZ", "BD", "Converter", "Non-Converter")) +  # Use `scale_fill_manual` for discrete data
    labs(
    title = "Proportions of CpG-Associated SCZ Cases by Region",
    x = "Region",
    y = "Proportion",
    fill = "Category"
    ) +
    theme_bw() +
    theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
    )
    print(stacked_bar_SCZ_CpG)



# all proportion for scz
    heatmap_data_SCZ <- results_SCZ %>%
    select(region, prop_SCZ_pure, prop_BD_pure, prop_Converter_pure, prop_Non_Converter_pure, prop_CpG_SCZ, prop_CpG_BD, prop_CpG_Converter, prop_CpG_Non_Converter) %>%
    pivot_longer(
    cols = starts_with("prop_"),
    names_to = "category",
    values_to = "proportion"
    )
    # Plot
    heatmap_SCZ <- ggplot(heatmap_data_SCZ, aes(x = region, y = category, fill = proportion)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(
    title = "Proportions of SCZ-Associated Cases by Region",
    x = "Region",
    y = "Category",
    fill = "Proportion"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(heatmap_SCZ)

















#
# krusko
    # Function to create boxplots with Dunn test results
    create_boxplot_with_dunn <- function(data, response, group, kruskal_result, dunn_test_result , name ) {
        p <- ggplot(data, aes_string(x = group, y = response, fill = group)) +
            geom_boxplot() +
            scale_fill_brewer(palette = "Set3") +
            labs(title = paste("Boxplot of", name, "by", group),
                x = group,
                y = response) +
            theme_minimal() +
            theme(legend.position = "none")
        
        if (!is.null(dunn_test_result)) {
            # Extract comparisons and p-values
            comparisons <- dunn_test_result$comparisons
            p_values <- dunn_test_result$altP.adjusted
            
            # Create a data frame for the annotations
            annotations <- data.frame(
                comparisons = comparisons,
                p_values = p_values
            )
            
            # Add significance stars to the annotations
            annotations <- annotations %>%
                mutate(stars = case_when(
                    p_values < 0.001 ~ "***",
                    p_values < 0.01 ~ "**",
                    p_values < 0.05 ~ "*",
                    TRUE ~ "ns"
                ))
            
            # Add annotations to the plot
            for (i in 1:nrow(annotations)) {
                comparison <- unlist(strsplit(annotations$comparisons[i], " - "))
                p <- p + geom_signif(
                    comparisons = list(comparison),
                    map_signif_level = TRUE,
                    annotations = annotations$stars[i],
                    y_position = max(data[[response]], na.rm = TRUE) + (i * 0.1) * max(data[[response]], na.rm = TRUE),
                    tip_length = 0.03
                )
            }
        }
        
        # Add Kruskal-Wallis result
        p <- p + annotate("text", x = 1, y = Inf, label = paste("K-W p =", formatC(kruskal_result$p.value, format = "e", digits = 2)), vjust = 2, hjust = -0.1)
        print(p)
        return(p)
    }


    kruskal_SCZ <-create_boxplot_with_dunn(data_tr, "tr_count", "group", kruskal_result$kruskal_result, kruskal_result$dunn_test_result , "TRs")
    kruskal_SCZ_normalized <-create_boxplot_with_dunn(data_long_normalized, "tr_per_individual", "group", kruskal_result_normalized$kruskal_result, kruskal_result_normalized$dunn_test_result , "TRs (normalized)")
    

    # do one pro region 
    for (region in unique(data_tr$region)) {
        # Get the results for the current region from the list
        region_result <- krusk_results_per_region[[region]]
        
        # Extract Kruskal-Wallis and Dunn results
        kruskal_res <- region_result$kruskal_result
        dunn_res <- region_result$dunn_test_result
        
        # Filter data for the current region
        region_data <- data_tr %>% filter(region == !!region)
        
        # Create the boxplot
        create_boxplot_with_dunn(
            region_data,
            "tr_count",
            "group",
            kruskal_res,
            dunn_res, region
        )
    }


    #FOR ALL REGIONS grid
    # Initialize an empty list for plots
    plots <- list()

    # Loop through each region and create plots
    for (region in unique(data_tr$region)) {
    # Filter data for the current region
    region_data <- data_tr %>% filter(region == !!region)
    
    # Get results for the current region
    region_result <- krusk_results_per_region[[region]]
    
    # Extract Kruskal-Wallis and Dunn results
    kruskal_res <- region_result$kruskal_result
    dunn_res <- region_result$dunn_test_result
    
    # Create the boxplot
    p <- create_boxplot_with_dunn(
        region_data,
        "tr_count",
        "group",
        kruskal_res,
        dunn_res,region
    )
    
    # Store the plot in the list
    plots[[region]] <- p
    }

    # Split the plots into chunks of 4 for each slide/page
    plot_chunks <- split(plots, ceiling(seq_along(plots) / 4))

    # Print each chunk of 4 plots on a new page
    for (chunk in plot_chunks) {
    grid.arrange(
        grobs = chunk,
        ncol = 2,
        top = "Boxplots with Dunn Test Annotations"
    )
    }

    #for all regions normalized
    # Initialize an empty list for plots
    plots_normalized <- list()
    # Loop through each region and create plots
    for (region in unique(data_tr$region)) {
    # Filter data for the current region
    region_data <- data_long_normalized %>% filter(region == !!region)
    # Get results for the current region
    region_result <- krusk_results_per_region_normalized[[region]]
    # Extract Kruskal-Wallis and Dunn results
    kruskal_res <- region_result$kruskal_result
    dunn_res <- region_result$dunn_test_result
    # Create the boxplot
    p <- create_boxplot_with_dunn(
        region_data,
        "tr_per_individual",
        "group",
        kruskal_res,
        dunn_res,region
    )
    # Store the plot in the list
    plots_normalized[[region]] <- p
    }
    # Split the plots into chunks of 4 for each slide/page
    plot_chunks_normalized <- split(plots_normalized, ceiling(seq_along(plots_normalized) / 4))
    # Print each chunk of 4 plots on a new page
    for (chunk in plot_chunks_normalized) {
    grid.arrange(
        grobs = chunk,
        ncol = 2,
        top = "Boxplots with Dunn Test Annotations (Normalized)"
    )
    }



#
# mann witheny 
# for region unnormalized 
    annotation_data <- wilcox_results_per_region %>%
        mutate(
            label = paste("FDR p =", format.pval(adjusted_p_value_fdr, digits = 2))
        )

    # 2. Create the plot
    plot_mann<-ggplot(data_tr_case, aes(x = group, y = tr_count)) +
        # Boxplots with jittered points
        geom_boxplot(aes(fill = group), outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
        
        # Split into facets by region
        facet_wrap(~ region, scales = "free_y") +
        
        # Add adjusted p-values as annotations
        geom_text(
            data = annotation_data,
            aes(x = 1.5, y = Inf, label = label),  # Center text between groups
            vjust = 1.5, size = 3, color = "red"
        ) +
        
        # Customize labels and theme
        labs(
            title = "Case vs. Control tr_count by Region",
            x = "Group",
            y = "Expansion amount (tr_count)",
            fill = "Group"
        ) +
        theme_bw() +
        theme(
            strip.background = element_blank(),
            strip.text = element_text(face = "bold")
        ) 
# for region normalized
    annotation_data <- wilcox_results_per_region_normalized %>%
        mutate(
            label = paste("FDR p =", format.pval(adjusted_p_value_fdr, digits = 2))
        )
        # 2. Create the plot
        plot_mann_normalized<-ggplot(data_long_normalized_case, aes(x = group, y = tr_per_individual)) +
        # Boxplots with jittered points
        geom_boxplot(aes(fill = group), outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.3, size = 1) +

        # Split into facets by region
        facet_wrap(~ region, scales = "free_y") +

        # Add adjusted p-values as annotations
        geom_text(
            data = annotation_data,
            aes(x = 1.5, y = Inf, label = label),  # Center
            vjust = 1.5, size = 3, color = "red"
            ) +
            # Customize labels and theme
            labs(
                title = "Case vs. Control tr_per_individual by Region (normalized)",
                x = "Group",
                y = "Expansion amount (tr_per_individual)",
                fill = "Group"
            ) +   
            theme_bw() +
            theme(
                strip.background = element_blank(),
                strip.text = element_text(face = "bold")
            )

# VISUALIZATION
    count_data <- results %>%
        select(region, total, case_pure, control_pure, detail_case, detail_control, total_unique_cases, total_unique_controls) %>%
        pivot_longer(cols = -region, names_to = "metric", values_to = "count")

    ggplot(count_data, aes(x = region, y = count, fill = metric)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "Counts of TRs and Associations by Region", x = "Region", y = "Count", fill = "Metric") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    association_data <- results %>%
        select(region, prop_cases_associated, prop_unique_cases_associated, prop_controls_associated, prop_unique_controls_associated) %>%
        pivot_longer(cols = -region, names_to = "metric", values_to = "proportion")

    ggplot(association_data, aes(x = region, y = proportion, fill = metric)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "Proportions of Cases and Controls Associated with TRs", x = "Region", y = "Proportion", fill = "Metric") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    association_data_SCZ <- results_SCZ %>%
        select(region, prop_SCZ_associated, prop_unique_SCZ_associated, prop_BD_associated, prop_unique_BD_associated, prop_Converter_associated, prop_unique_Converter_associated, prop_Non_Converter_associated, prop_unique_Non_Converter_associated) %>%
        pivot_longer(cols = -region, names_to = "metric", values_to = "proportion")

    ggplot(association_data_SCZ, aes(x = region, y = proportion, fill = metric)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "Proportions of SCZ, BD, Converter, and Non-Converter Associated with TRs", x = "Region", y = "Proportion", fill = "Metric") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))


    count_data_SCZ <- results_SCZ %>%
        select(region, total, SCZ_pure, BD_pure, Converter_pure, Non_Converter_pure, detail_SCZ, detail_BD, detail_converter, detail_non_converter, total_unique_SCZ, total_unique_BD, total_unique_Converter, total_unique_Non_Converter) %>%
        pivot_longer(cols = -region, names_to = "metric", values_to = "count")

    ggplot(count_data_SCZ, aes(x = region, y = count, fill = metric)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "Counts of TRs and Associations by Region (SCZ, BD, Converter, Non-Converter)", x = "Region", y = "Count", fill = "Metric") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Summary tables
    summary_table <- results %>%
        select(region, prop_pure_case, prop_case_mix, prop_pure_control, prop_control_mix, prop_cases_associated, prop_unique_cases_associated, prop_controls_associated, prop_unique_controls_associated) %>%
        rename(
            "Region" = region,
            "Proportion in Pure Cases" = prop_pure_case,
            "Proportion in Cases and Mixed" = prop_case_mix,
            "Proportion in Pure Controls" = prop_pure_control,
            "Proportion in Controls and Mixed" = prop_control_mix,
            "Proportion of Cases Associated" = prop_cases_associated,
            "Proportion of Unique Cases Associated" = prop_unique_cases_associated,
            "Proportion of Controls Associated" = prop_controls_associated,
            "Proportion of Unique Controls Associated" = prop_unique_controls_associated
        )

    summary_table %>%
        kable("html", caption = "Summary of TR Proportions and Associations") %>%
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
        column_spec(1, bold = TRUE) %>%
        add_header_above(c(" " = 1, "Cases" = 3, "Controls" = 3, "Associations" = 2))

    summary_table_SCZ <- results_SCZ %>%
        select(region, prop_SCZ_pure, prop_SCZ_mix, prop_BD_pure, prop_BD_mix, prop_Converter_pure, prop_Converter_mix, prop_Non_Converter_pure, prop_Non_Converter_mix, prop_SCZ_associated, prop_unique_SCZ_associated, prop_BD_associated, prop_unique_BD_associated, prop_Converter_associated, prop_unique_Converter_associated, prop_Non_Converter_associated, prop_unique_Non_Converter_associated) %>%
        rename(
            "Region" = region,
            "Proportion in Pure SCZ" = prop_SCZ_pure,
            "Proportion in SCZ and Mixed" = prop_SCZ_mix,
            "Proportion in Pure BD" = prop_BD_pure,
            "Proportion in BD and Mixed" = prop_BD_mix,
            "Proportion in Pure Converter" = prop_Converter_pure,
            "Proportion in Converter and Mixed" = prop_Converter_mix,
            "Proportion in Pure Non-Converter" = prop_Non_Converter_pure,
            "Proportion in Non-Converter and Mixed" = prop_Non_Converter_mix,
            "Proportion of SCZ Associated" = prop_SCZ_associated,
            "Proportion of Unique SCZ Associated" = prop_unique_SCZ_associated,
            "Proportion of BD Associated" = prop_BD_associated,
            "Proportion of Unique BD Associated" = prop_unique_BD_associated,
            "Proportion of Converter Associated" = prop_Converter_associated,
            "Proportion of Unique Converter Associated" = prop_unique_Converter_associated,
            "Proportion of Non-Converter Associated" = prop_Non_Converter_associated,
            "Proportion of Unique Non-Converter Associated" = prop_unique_Non_Converter_associated
        )

    summary_table_SCZ %>%
        kable("html", caption = "Summary of TR Proportions and Associations (SCZ, BD, Converter, Non-Converter)") %>%
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
        column_spec(1, bold = TRUE) %>%
        add_header_above(c(" " = 1, "SCZ" = 2, "BD" = 2, "Converter" = 2, "Non-Converter" = 2, "Associations" = 8))




# create tabel for kruskal for region
    process_kruskal_results <- function(kruskal_list, normalized = FALSE) {
        # Process Kruskal-Wallis results
        kruskal_df <- map_dfr(kruskal_list, function(region_results) {
            tibble(
                region = region_results$kruskal_result$data.name,
                statistic = region_results$kruskal_result$statistic,
                p.value = region_results$kruskal_result$p.value,
                method = "Kruskal-Wallis",
                normalized = normalized
            )
        }, .id = "region")
        
        # Process Dunn's test results if present
        dunn_df <- map_dfr(kruskal_list, function(region_results) {
            if (!is.null(region_results$dunn_test_result)) {
                tibble(
                    region = region_results$kruskal_result$data.name,
                    comparison = region_results$dunn_test_result$comparisons,
                    Z = region_results$dunn_test_result$Z,
                    p.value = region_results$dunn_test_result$altP,
                    p.adj = region_results$dunn_test_result$altP.adjusted,
                    method = "Dunn's test",
                    normalized = normalized
                )
            }
        }, .id = "region")
        
        list(kruskal = kruskal_df, dunn = dunn_df)
    }

    # Process non-normalized results
    kruskal_non_norm <- process_kruskal_results(krusk_results_per_region, normalized = FALSE)
    kruskal_norm <- process_kruskal_results(krusk_results_per_region_normalized, normalized = TRUE)


    # Combine all results
    combined_kruskal_region <- bind_rows(
        kruskal_non_norm$kruskal,
        kruskal_norm$kruskal,
        kruskal_non_norm$dunn,
        kruskal_norm$dunn
    ) 
# create tabel for kruskale no region 

    # Function to process Kruskal-Wallis and Dunn's test results
    process_kruskal_dunn_results <- function(kruskal_dunn_result, normalized = FALSE) {
    # Extract Kruskal-Wallis results
    kruskal_df <- tibble(
        region = "Overall",  # Since this is not region-specific, use "Overall"
        statistic = kruskal_dunn_result$kruskal_result$statistic,
        p.value = kruskal_dunn_result$kruskal_result$p.value,
        method = "Kruskal-Wallis",
        normalized = normalized
    )
    
    # Extract Dunn's test results if present
    dunn_df <- NULL
    if (!is.null(kruskal_dunn_result$dunn_test_result)) {
        dunn_df <- tibble(
        region = "Overall",
        comparison = kruskal_dunn_result$dunn_test_result$comparisons,
        Z = kruskal_dunn_result$dunn_test_result$Z,
        p.value = kruskal_dunn_result$dunn_test_result$altP,
        p.adj = kruskal_dunn_result$dunn_test_result$altP.adjusted,
        method = "Dunn's test",
        normalized = normalized
        )
    }
    
    list(kruskal = kruskal_df, dunn = dunn_df)
    }

    # Process non-normalized results
    kruskal_non_norm <- process_kruskal_dunn_results(kruskal_result, normalized = FALSE)

    # Process normalized results
    kruskal_norm <- process_kruskal_dunn_results(kruskal_result_normalized, normalized = TRUE)

    # Combine all results
    combined_kruskal <- bind_rows(
    kruskal_non_norm$kruskal,
    kruskal_norm$kruskal,
    kruskal_non_norm$dunn,
    kruskal_norm$dunn
    )

    # Print the combined table
    print(combined_kruskal)


#
# save plots for Psycho vs Case
# Save Psycho vs Control Results
    setwd("/home/rachele/EHDN_DBSCAN_correct/Result/results_dbscan_with_QC_NOUHR/AFTER_DBSCAN/psycho_vs_control")

    # Save Plots
    ggsave("manhattan_plot.png", manhattan_plot, width=10, height=8)
    ggsave("combined_volcano.png", combined_volcano, width=12, height=6)
    ggsave("vulcanoplot_tr.png", volcano_tr, width=10, height=8)
    ggsave("dot_plot.png", dot_plot, width=10, height=8)
    ggsave("stacked_bar_pure.png", stacked_bar_pure, width=12, height=8)
    ggsave("stacked_bar_cpg.png", stacked_bar_cpg, width=12, height=8)
    ggsave("heatmap.png", heatmap, width=15, height=13)
    ggsave("mann_whitney_plots.png", plot_mann, width=14, height=10)
    ggsave("mann_whitney_normalized.png", plot_mann_normalized, width=14, height=10)

    # Save Tables
    write.csv(results, "case_control_summary.csv", row.names=FALSE)
    write.csv(counts_case_control, "counts_case_control.csv", row.names=FALSE)
    write.csv(fisher_results_tr, "fisher_test_TR.csv", row.names=FALSE)
    write.csv(fisher_results_cpg, "fisher_test_CpG.csv", row.names=FALSE)
    write.csv(fisher_results, "fisher_test_individuals.csv", row.names=FALSE)
    write.csv(wilcox_results_per_region, "mann_whitney_results.csv", row.names=FALSE)
    write.csv(wilcox_results_per_region_normalized, "mann_whitney_normalized_results.csv", row.names=FALSE)
    write.csv(wilcoxon_results_df, "wilcoxon_results_general.csv", row.names=FALSE)
    write.csv(wilcoxon_results_normalized_df, "wilcoxon_results_general_normalized.csv", row.names=FALSE)
    write.csv(summary_table, "summary_table.csv", row.names=FALSE)


# Save SCZ vs Control Results

setwd("/home/rachele/EHDN_DBSCAN_correct/Result/results_dbscan_with_QC_NOUHR/AFTER_DBSCAN/SCZvsControl")
#save plots
    ggsave("dot_plot_SCZ.png", dot_plot_SCZ, width=10, height=8)
    ggsave("vulcano_plot_SCZ.png", vulcano_plot_SCZ, width=12, height=6)
    ggsave("heatmap_SCZ_pure.png", pure_plot, width=10, height=8)
    ggsave("heatmap_SCZ_CpG.png", cpg_plot, width=10, height=8)
    ggsave("heatmap_individual.png", individual_plot, width=10, height=8)
    ggsave("combined_heatmap.png", heatmap_combined, width=19, height=13)
    ggsave("lollipop_SCZ.png", lolli_SCZ, width=10, height=8)
    ggsave("dot_plot_size_effect_SCZ.png", dot_plot_size_effect_SCZ, width=10, height=8)
    ggsave("stacked_bar_SCZ.png", stacked_bar_SCZ, width=12, height=8)
    ggsave("stacked_bar_SCZ_CpG.png", stacked_bar_SCZ_CpG, width=12, height=8)
    ggsave("heatmap_SCZ_all.png", heatmap_SCZ, width=10, height=8)
    ggsave("kruskal_SCZ.png", kruskal_SCZ, width=10, height=8)
    ggsave("kruskal_SCZ_normalized.png", kruskal_SCZ_normalized, width=10, height=8)


    # Save grid of plots (4 per page)
    # Split the plots into chunks of 4 for each page
    plot_chunks <- split(plots, ceiling(seq_along(plots) / 4))
    plot_chunks_normalized <- split(plots_normalized, ceiling(seq_along(plots_normalized) / 4))

    # Save raw data plots
    for (i in seq_along(plot_chunks)) {
    grid_plot <- grid.arrange(
        grobs = plot_chunks[[i]],
        ncol = 2,
        top = paste("Boxplots with Dunn Test Annotations (Raw Data) - Page", i)
    )
    ggsave(paste0("raw_data_grid_page_", i, ".png"), grid_plot, width = 10, height = 8)
    }

    # Save normalized data plots
    for (i in seq_along(plot_chunks_normalized)) {
    grid_plot <- grid.arrange(
        grobs = plot_chunks_normalized[[i]],
        ncol = 2,
        top = paste("Boxplots with Dunn Test Annotations (Normalized Data) - Page", i)
    )
    ggsave(paste0("normalized_data_grid_page_", i, ".png"), grid_plot, width = 10, height = 8)
    }

#Save tabels for these counts_SCZ , results_SCZ,summary_table_pure,summary_table_CpG,summary_table_FISHER_SCZ,summary_table_SCZ,combined_kruskal
write.csv(counts_SCZ, "counts_SCZ.csv", row.names=FALSE)
write.csv(results_SCZ, "results_SCZ.csv", row.names=FALSE)
write.csv(summary_table_pure, "summary_table_pure.csv", row.names=FALSE)
write.csv(summary_table_CpG, "summary_table_CpG.csv", row.names=FALSE)
write.csv(summary_table_FISHER_SCZ, "summary_table_FISHER_SCZ.csv", row.names=FALSE)
write.csv(summary_table_SCZ, "summary_table_SCZ.csv", row.names=FALSE)
write.csv(combined_kruskal, "combined_kruskal.csv", row.names=FALSE)
write.csv(combined_kruskal_region, "combined_kruskal_region.csv", row.names=FALSE)





#
# Bad 

    proportion_data <- results %>%
        select(region, prop_pure_case, prop_case_mix, prop_pure_control, prop_control_mix) %>%
        pivot_longer(cols = -region, names_to = "metric", values_to = "proportion")

    ggplot(proportion_data, aes(x = region, y = proportion, fill = metric)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "Proportions of TRs in Cases and Controls", x = "Region", y = "Proportion", fill = "Metric") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    proportion_data_SCZ <- results_SCZ %>%
        select(region, prop_SCZ_pure, prop_SCZ_mix, prop_BD_pure, prop_BD_mix, prop_Converter_pure, prop_Converter_mix, prop_Non_Converter_pure, prop_Non_Converter_mix) %>%
        pivot_longer(cols = -region, names_to = "metric", values_to = "proportion")

    ggplot(proportion_data_SCZ, aes(x = region, y = proportion, fill = metric)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "Proportions of TRs in SCZ, BD, Converter, and Non-Converter", x = "Region", y = "Proportion", fill = "Metric") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    # VISUALIZATION
        # Plot p-values by region
        ggplot(fisher_results_tr, aes(x = region, y = p.value, color = p.value < 0.05)) +
            geom_point(size = 3) +
            geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
            theme_minimal() +
            labs(title = "P-values by Region (TR-Based Analysis)", x = "Region", y = "P-value") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))

        # Plot p-values by region for CpG-related TRs
        ggplot(fisher_results_cpg, aes(x = region, y = p.value, color = p.value < 0.05)) +
            geom_point(size = 3) +
            geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
            theme_minimal() +
            labs(title = "P-values by Region (CpG-Related TRs)", x = "Region", y = "P-value") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))
# ELIMINATE 
# Kruskal 
    # Reshape the data for Kruskal-Wallis test
    long_data <- results_SCZ %>%
        select(region, SCZ_pure, BD_pure, Converter_pure, Non_Converter_pure) %>%
        pivot_longer(
            cols = c(SCZ_pure, BD_pure, Converter_pure, Non_Converter_pure),
            names_to = "group",
            values_to = "count"
        )

    kruskal_test_result <- kruskal.test(count ~ group, data = long_data)

    dunn_test_result <- dunnTest(count ~ group, data = long_data, method = "bh")

    ggplot(long_data, aes(x = group, y = count, fill = group)) +
        geom_boxplot() +
        theme_minimal() +
        labs(title = "Distribution of TR Counts by Group", x = "Group", y = "TR Count")



    # Reshape the data for CpG-related TRs
    long_data_cpg <- results_SCZ %>%
        select(region, CpG_SCZ, CpG_BD, CpG_Converter, CpG_Non_Converter) %>%
        pivot_longer(
            cols = c(CpG_SCZ, CpG_BD, CpG_Converter, CpG_Non_Converter),
            names_to = "group",
            values_to = "count"
        )

    # Perform Kruskal-Wallis test for CpG-related TRs
    kruskal_test_result_cpg <- kruskal.test(count ~ group, data = long_data_cpg)

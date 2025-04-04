#statistics 3 groups 


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


    ehdn_private <- ehdn_results_annotated %>%filter(count == 1 )



#prepocess SCZ 
counts_SCZ <- ehdn_private %>% 
  group_by(region) %>%
  summarise(
    total = n(),
    FEP_pure = sum(outlier_label_detail %in% c("FEP-SCZ", "FEP-BD")),
    Converter_pure = sum(outlier_label_detail == "Converter"),
    Non_Converter_pure = sum(outlier_label_detail == "Non_Converter"),
    detail_FEP = sum(count_SCZ + count_BD),  # Combined FEP counts
    detail_converter = sum(count_converter),
    detail_non_converter = sum(count_non_converter),
    CpG_FEP = sum(CpG[outlier_label_detail %in% c("FEP-SCZ", "FEP-BD")]),
    CpG_Converter = sum(CpG[outlier_label_detail == "Converter"]),
    CpG_Non_Converter = sum(CpG[outlier_label_detail == "Non_Converter"]),
    CpG = sum(CpG)
  )



# Read manifest file and expand 'outliers' into separate rows
manifest <- read_excel("~/Downloads/Liste_samples_WGS_final.xlsx")

# Combine FEP-SCZ and FEP-BD counts for unique samples
unique_case_control_counts <- ehdn_private %>%
  separate_rows(outliers, sep = ";") %>%
  left_join(manifest, by = c("outliers" = "Sequencing_number")) %>%
  mutate(
    Status = case_when(
      Status %in% c("FEP-SCZ", "FEP-BD") ~ "FEP",  # Combine FEP groups
      TRUE ~ Status  # Keep other statuses as-is
    )
  ) %>%
  group_by(region, Status) %>%
  summarise(unique_samples = n_distinct(outliers), .groups = "drop") %>%
  pivot_wider(
    names_from = Status, 
    values_from = unique_samples, 
    values_fill = list(unique_samples = 0)
  ) %>%
  rename(
    total_unique_FEP = FEP,
    total_unique_Converter = Converter,
    total_unique_Non_Converter = Non_Converter
  )

# Merge with main summary counts
counts_SCZ <- counts_SCZ %>%
  left_join(unique_case_control_counts, by = "region")

# Calculate proportions and counts
total_FEP <- 255  # Sum of total_SCZ (222) + total_BD (33)
total_converter <- 47
total_non_converter <- 75

# Calculate proportions and counts for SCZ
results_SCZ <- counts_SCZ %>%
  mutate(
    prop_FEP_pure = FEP_pure / total,
    prop_Converter_pure = Converter_pure / total,
    prop_Non_Converter_pure = Non_Converter_pure / total,
    prop_FEP_associated = detail_FEP / total_FEP,
    prop_unique_FEP_associated = total_unique_FEP / total_FEP,
    prop_Converter_associated = detail_converter / total_converter,
    prop_unique_Converter_associated = total_unique_Converter / total_converter,
    prop_Non_Converter_associated = detail_non_converter / total_non_converter,
    prop_unique_Non_Converter_associated = total_unique_Non_Converter / total_non_converter,
    prop_CpG_FEP = CpG_FEP / total,
    prop_CpG_Converter = CpG_Converter / total,
    prop_CpG_Non_Converter = CpG_Non_Converter / total,
    prop_CpG = CpG / total
  )


# FISHER _INDIVIDUALS
# Update Fisher tests for individuals with combined FEP group
pairwise_data <- counts_SCZ %>%
  select(region, FEP_pure, Converter_pure, Non_Converter_pure) %>%
  mutate(
    FEP_non_outlier = total_FEP - FEP_pure,
    Converter_non_outlier = total_converter - Converter_pure,
    Non_Converter_non_outlier = total_non_converter - Non_Converter_pure
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

# Define the new groups to compare
groups_to_compare <- c("FEP", "Converter", "Non_Converter")

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




#fisher individuals cpg 
pairwise_data <- counts_SCZ %>%
  select(region, CpG_FEP, CpG_Converter,CpG_Non_Converter) %>%
  mutate(
    FEP_non_CpG = total_FEP - CpG_FEP,
    Converter_non_CpG = total_converter - CpG_Converter,
    Non_Converter_non_CpG = total_non_converter - CpG_Non_Converter
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
                            group1_CpG = paste0( "CpG_",group1),
                            group2_CpG = paste0( "CpG_",group2),
                            group1_non_CpG = paste0(group1, "_non_CpG"),
                            group2_non_CpG = paste0(group2, "_non_CpG")
                        )
                    
                    # Perform Fisher's exact test for each region
                    fisher_results <- contingency_table %>%
                        group_by(region) %>%
                        group_modify(~ broom::tidy(fisher.test(matrix(c(
                            .x$group1_CpG, .x$group2_CpG,  # Group 1 and Group 2 counts
                            .x$group1_non_CpG, .x$group2_non_CpG  # Non-outlier counts
                        ), nrow = 2))))
                        
                    # Store the results
                    results[[paste(group1, "vs", group2)]] <- fisher_results
                }
            }
            
            return(results)
        }

        # Define the new groups to compare
        groups_to_compare <- c("FEP", "Converter", "Non_Converter")


            # Perform pairwise Fisher's exact tests
            pairwise_results_CpG <- perform_pairwise_fisher_individuals(pairwise_data, groups_to_compare)

            summary_table_FISHER_SCZ_CpG <- do.call(rbind, lapply(names(pairwise_results_CpG), function(comparison) {
                pairwise_results_CpG[[comparison]] %>%
                    mutate(Comparison = comparison)
            }))

            summary_table_FISHER_SCZ_CpG <- summary_table_FISHER_SCZ_CpG %>%
                mutate(
                    adjusted_p_value_fdr = p.adjust(p.value, method = "fdr"),
                    adjusted_p_value_bonferroni = p.adjust(p.value, method = "bonferroni")
                )



#kruskal 
#for normal TREs
    data_tr <- ehdn_private %>%
    select(region, count_SCZ, count_BD, count_converter, count_non_converter) %>%
    pivot_longer(
        cols = c(count_SCZ, count_BD, count_converter, count_non_converter),
        names_to = "group",
        values_to = "tr_count"
    ) %>%
    mutate(
        group = case_when(
        group %in% c("count_SCZ", "count_BD") ~ "FEP",  # Combine SCZ and BD
        group == "count_converter" ~ "Converter",
        group == "count_non_converter" ~ "Non_Converter",
        TRUE ~ "Other"
        )
    ) %>%
    group_by(region, group) %>%
    summarise(tr_count = sum(tr_count), .groups = "drop")  # Sum counts for combined FEP group

# For CpG analysis with combined FEP group
    data_tr_cpg <- ehdn_private %>%
    select(region, outlier_label_detail, CpG) %>%
    filter(!is.na(outlier_label_detail)) %>%
    mutate(
        group = case_when(
        outlier_label_detail %in% c("FEP-SCZ", "FEP-BD") ~ "FEP",  # Combine SCZ and BD
        outlier_label_detail == "Converter" ~ "Converter",
        outlier_label_detail == "Non_Converter" ~ "Non_Converter",
        TRUE ~ "Other"
        )
    ) %>%
    select(region, group, tr_count = CpG)

#automatices
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



    #new with more safeguards 
    perform_kruskal_dunn_test_region <- function(data, response_col, group_col) {
  results <- list()
  
  for (region in unique(data$region)) {
    region_data <- data %>% filter(region == !!region)
    
    # Check if there are at least 2 groups in the region
    unique_groups <- unique(region_data[[group_col]])
    if (length(unique_groups) < 2) {
      warning(paste("Skipping region", region, ": Only one group (", unique_groups, ") present."))
      next  # Skip to the next region
    }
    
    kruskal_formula <- as.formula(paste(response_col, "~", group_col))
    kruskal_result <- kruskal.test(kruskal_formula, data = region_data)
    print(kruskal_result)
    
    dunn_test_result <- NULL
    if (!is.na(kruskal_result$p.value) && kruskal_result$p.value <= 0.05) {
      dunn_test_result <- dunn.test(region_data[[response_col]], region_data[[group_col]], method = "bh", altp = TRUE)
      print(dunn_test_result)
    } else {
      print(paste("Skipping Dunn's test for region", region, "due to invalid or non-significant Kruskal-Wallis p-value"))
    }
    
    results[[region]] <- list(kruskal_result = kruskal_result, dunn_test_result = dunn_test_result)
  }
  
  return(results)
}


#execution
    # Non-normalized data
    kruskal_result <- perform_kruskal_dunn_test(data_tr, "tr_count", "group")
    krusk_results_per_region <- perform_kruskal_dunn_test_region(data_tr, "tr_count", "group")


    #for cpg
    kruskal_result_cpg <- perform_kruskal_dunn_test(data_tr_cpg, "tr_count", "group")
    krusk_results_per_region_cpg <- perform_kruskal_dunn_test_region(data_tr_cpg, "tr_count", "group")





# FOR SCZ 
# FISHER
# dot_plot
    dot_plot_SCZ <- ggplot(summary_table_FISHER_SCZ, aes(x = -log10(adjusted_p_value_fdr), y = region)) +
        geom_point(aes(color = Comparison , shape = adjusted_p_value_fdr < 0.05), size = 3) +
        geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
        
        scale_shape_manual(values = c(16, 17)) +  # Shapes for non-sig vs. sig
        labs(
            title = "Significance by Region and comparison",
            x = "-log10(FDR-adjusted p-value)",
            y = "Region",
            color = "comparison",
            shape = "Significant (FDR < 0.05)"
        ) +
        theme_bw() +
        theme(axis.text.y = element_text(size = 8))

    print(dot_plot_SCZ)

# Volcano plot for TR outliers

    summary_table_FISHER_SCZ  <- summary_table_FISHER_SCZ %>%
    mutate(
        log_OR = log2(estimate),
        log_p = -log10(adjusted_p_value_fdr),
        Significance = ifelse(adjusted_p_value_fdr < 0.05, "FDR < 0.05", "NS")
    )

    vulcano_plot_SCZ <- ggplot(summary_table_FISHER_SCZ, aes(x = log_OR, y = log_p, color = Significance)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.8, size = 3) +
    scale_color_manual(values = c("FDR < 0.05" = "red", "NS" = "grey60")) +
    theme_bw() +
    labs(
        x = "Log2(Odds Ratio)", 
        y = "-Log10(FDR-adjusted p)",
        title = "Volcano Plots of Pairwise Comparisons"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Heatmap of log odds ratios + p value

    heatmap_fisher_indiv <- ggplot(summary_table_FISHER_SCZ, aes(x = Comparison, y = region, fill = log_OR)) +
        geom_tile(color = "white") +
        geom_text(aes(label = ifelse(Significance == "FDR < 0.05", "*", "")), 
                vjust = 0.8, size = 6) +
        scale_fill_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = 0, name = "Log2(OR)"
        ) +
        theme_minimal() +
        theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")
        )
  


# lollipop plot for effect size 
    lolli_SCZ <-ggplot(summary_table_FISHER_SCZ, aes(x = log_OR, y = reorder(region, log_OR))) +
    geom_segment(aes(xend = 0, yend = region), color = "grey50") +
    geom_point(aes(
        color = ifelse(adjusted_p_value_fdr < 0.05, "Sig", "Non-Sig"),
        size = 3
    )) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("Sig" = "red", "Non-Sig" = "grey60")) +
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

# dot_plot  for size effect and Direction
    plot_data <- summary_table_FISHER_SCZ %>%
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
    select(region, prop_FEP_pure, prop_Converter_pure, prop_Non_Converter_pure) %>%
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
    scale_fill_manual(values = scz_colors, labels = c("FEP", "Converter", "Non-Converter")) +  # Use `scale_fill_manual` for discrete data
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
    select(region, prop_CpG_FEP, prop_CpG_Converter, prop_CpG_Non_Converter) %>%
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
    scale_fill_manual(values = scz_colors, labels = c("FEP", "Converter", "Non-Converter")) +  # Use `scale_fill_manual` for discrete data
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
    select(region, prop_FEP_pure, prop_Converter_pure, prop_Non_Converter_pure, prop_CpG_FEP, prop_CpG_Converter, prop_CpG_Non_Converter) %>%
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

# general visualization 
    association_data_SCZ <- results_SCZ %>%
        select(region, prop_FEP_associated, prop_unique_FEP_associated, prop_Converter_associated, prop_unique_Converter_associated, prop_Non_Converter_associated, prop_unique_Non_Converter_associated) %>%
        pivot_longer(cols = -region, names_to = "metric", values_to = "proportion")

    proportion_plot_SCZ <- ggplot(association_data_SCZ, aes(x = region, y = proportion, fill = metric)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "Proportions of FEP, Converter, and Non-Converter Associated with TRs", x = "Region", y = "Proportion", fill = "Metric") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))


    count_data_SCZ <- results_SCZ %>%
        select(region, total, FEP_pure,  Converter_pure, Non_Converter_pure, total_unique_FEP, total_unique_Converter, total_unique_Non_Converter) %>%
        pivot_longer(cols = -region, names_to = "metric", values_to = "count")

    count_plot_SCZ <- ggplot(count_data_SCZ, aes(x = region, y = count, fill = metric)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "Counts of TRs and Associations by Region (FEP, Converter, Non-Converter)", x = "Region", y = "Count", fill = "Metric") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))




#summary tabel 

    
    summary_table_SCZ <- results_SCZ %>%
        select(region, prop_FEP_pure,  prop_Converter_pure, prop_Non_Converter_pure, prop_FEP_associated, prop_unique_FEP_associated, prop_Converter_associated, prop_unique_Converter_associated, prop_Non_Converter_associated, prop_unique_Non_Converter_associated) %>%
        rename(
            "Region" = region,
            "Proportion in Pure FEP" = prop_FEP_pure,
            "Proportion in Pure Converter" = prop_Converter_pure,
            "Proportion in Pure Non-Converter" = prop_Non_Converter_pure,
            "Proportion of FEP Associated" = prop_FEP_associated,
            "Proportion of Unique FEP Associated" = prop_unique_FEP_associated,
            "Proportion of Converter Associated" = prop_Converter_associated,
            "Proportion of Unique Converter Associated" = prop_unique_Converter_associated,
            "Proportion of Non-Converter Associated" = prop_Non_Converter_associated,
            "Proportion of Unique Non-Converter Associated" = prop_unique_Non_Converter_associated
        )

    summary_table_SCZ %>%
        kable("html", caption = "Summary of TR Proportions and Associations (FEP, Converter, Non-Converter)") %>%
        kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
        column_spec(1, bold = TRUE) %>%
        add_header_above(c(" " = 1, "FEP" = 1,  "Converter" = 1, "Non-Converter" = 1, "Associations" = 6))    




# kruskal tabel 

  # tabel for kruskal for region
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


    # Combine all results
    combined_kruskal_region <- bind_rows(
        kruskal_non_norm$kruskal,
        kruskal_non_norm$dunn)


#  tabel for kruskale no region 

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

    # Combine all results
    combined_kruskal <- bind_rows(
    kruskal_non_norm$kruskal,
    kruskal_non_norm$dunn
    )

    # Print the combined table
    print(combined_kruskal)


#save 

setwd("/home/rachele/EHDN_DBSCAN_correct/Result/results_dbscan_with_QC_NOUHR/AFTER_DBSCAN/private/SCZvsControl/group3")
#save plots
    ggsave("dot_plot_SCZ.png", dot_plot_SCZ, width=10, height=8)
    ggsave("vulcano_plot_SCZ.png", vulcano_plot_SCZ, width=12, height=6)
    ggsave("heatmap_fisher_indiv.png", heatmap_fisher_indiv, width=10, height=8)
    ggsave("lollipop_SCZ.png", lolli_SCZ, width=10, height=8)
    ggsave("dot_plot_size_effect_SCZ.png", dot_plot_size_effect_SCZ, width=10, height=8)
    ggsave("stacked_bar_SCZ.png", stacked_bar_SCZ, width=12, height=8)
    ggsave("stacked_bar_SCZ_CpG.png", stacked_bar_SCZ_CpG, width=12, height=8)
    ggsave("heatmap_SCZ_all.png", heatmap_SCZ, width=10, height=8)
    ggsave("kruskal_SCZ.png", kruskal_SCZ, width=10, height=8)
    ggsave("count_plot_SCZ.png", count_plot_SCZ , width=14, height=10)
    ggsave("proportion_plot_SCZ.png", proportion_plot_SCZ, width=14, height=10)


    # Save grid of plots (4 per page)
    # Split the plots into chunks of 4 for each page
    plot_chunks <- split(plots, ceiling(seq_along(plots) / 4))

    # Save raw data plots
    for (i in seq_along(plot_chunks)) {
    grid_plot <- grid.arrange(
        grobs = plot_chunks[[i]],
        ncol = 2,
        top = paste("Boxplots with Dunn Test Annotations (Raw Data) - Page", i)
    )
    ggsave(paste0("raw_data_grid_page_", i, ".png"), grid_plot, width = 10, height = 8)
    }


#Save tabels for these counts_SCZ , results_SCZ,summary_table_pure,summary_table_CpG,summary_table_FISHER_SCZ,summary_table_SCZ,combined_kruskal
write.csv(counts_SCZ, "counts_SCZ.csv", row.names=FALSE)
write.csv(results_SCZ, "results_SCZ.csv", row.names=FALSE)
write.csv(summary_table_FISHER_SCZ_CpG, "summary_table_CpG.csv", row.names=FALSE)
write.csv(summary_table_FISHER_SCZ, "summary_table_FISHER_SCZ.csv", row.names=FALSE)
write.csv(summary_table_SCZ, "summary_table_SCZ.csv", row.names=FALSE)
write.csv(combined_kruskal, "combined_kruskal.csv", row.names=FALSE)
write.csv(combined_kruskal_region, "combined_kruskal_region.csv", row.names=FALSE)



# for only signifcant plots 
# kruskal 

create_boxplot_with_dunn <- function(data, response, group, kruskal_result, dunn_test_result, name) {
  # Only proceed if Kruskal-Wallis p-value is significant
  if (!is.null(kruskal_result$p.value) && kruskal_result$p.value < 0.05) {
    p <- ggplot(data, aes_string(x = group, y = response, fill = group)) +
      geom_boxplot() +
      scale_fill_brewer(palette = "Set3") +
      labs(title = paste("Boxplot of", name, "by", group),
           x = group,
           y = response) +
      theme_minimal() +
      theme(legend.position = "none")
    
    # Add Dunn's test annotations if available
    if (!is.null(dunn_test_result)) {
      comparisons <- dunn_test_result$comparisons
      p_values <- dunn_test_result$altP.adjusted
      
      annotations <- data.frame(
        comparisons = comparisons,
        p_values = p_values
      ) %>%
        mutate(stars = case_when(
          p_values < 0.001 ~ "***",
          p_values < 0.01 ~ "**",
          p_values < 0.05 ~ "*",
          TRUE ~ "ns"
        ))
      
      # Add only significant comparisons (p < 0.05)
      sig_annotations <- annotations %>% filter(p_values < 0.05)
      if (nrow(sig_annotations) > 0) {
        for (i in 1:nrow(sig_annotations)) {
          comparison <- unlist(strsplit(sig_annotations$comparisons[i], " - "))
          p <- p + geom_signif(
            comparisons = list(comparison),
            annotations = sig_annotations$stars[i],
            y_position = max(data[[response]], na.rm = TRUE) + (i * 0.1) * max(data[[response]], na.rm = TRUE),
            tip_length = 0.03
          )
        }
      }
    }
    
    # Add Kruskal-Wallis p-value
    p <- p + annotate("text", x = 1, y = Inf, 
                      label = paste("K-W p =", formatC(kruskal_result$p.value, format = "e", digits = 2)), 
                      vjust = 2, hjust = -0.1)
    return(p)
  } else {
    return(NULL)  # Return NULL if not significant
  }
}

# Initialize list for significant plots
sig_plots <- list()

for (region in unique(data_tr$region)) {
  region_data <- data_tr %>% filter(region == !!region)
  region_result <- krusk_results_per_region[[region]]
  p <- create_boxplot_with_dunn(region_data, "tr_count", "group", 
                               region_result$kruskal_result, 
                               region_result$dunn_test_result, 
                               region)
  if (!is.null(p)) {
    sig_plots[[region]] <- p
  }
}

# Plot significant regions in a grid (4 per page)
if (length(sig_plots) > 0) {
  plot_chunks <- split(sig_plots, ceiling(seq_along(sig_plots) / 4))
  for (chunk in plot_chunks) {
    grid.arrange(grobs = chunk, ncol = 2, top = "Significant Regions (Kruskal-Wallis p < 0.05)")
  }
} else {
  message("No significant regions found for Kruskal-Wallis.")
}


    # Save raw data plots
    for (i in seq_along(plot_chunks)) {
    grid_plot <- grid.arrange(
        grobs = plot_chunks[[i]],
        ncol = 2,
        top = paste("Boxplots with Dunn Test Annotations (Raw Data) - Page", i)
    )
    ggsave(paste0("Significant_raw_data_grid_page_", i, ".png"), grid_plot, width = 10, height = 8)
    }



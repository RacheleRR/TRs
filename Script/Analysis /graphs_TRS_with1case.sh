#!/bin/bash/
# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Read the reordered and annotated TSV files
ehdn_results_reordered <- read.delim("/home/rachele/ehdn_DBSCAN_reorder.tsv", stringsAsFactors = FALSE)
ehdn_results_annotated <- read.delim("/home/rachele/ehdn_DBSCAN_annotated.tsv", stringsAsFactors = FALSE)

# Merge reordered and annotated results
ehdn_results_annotated <- ehdn_results_annotated %>%
    left_join(ehdn_results_reordered %>% select(repeatID, ref, outlier_label, outlier_label2), by = "repeatID", suffix = c("", ".reordered")) %>%
    mutate(
        ref = coalesce(ref.reordered, ref),
        outlier_label = coalesce(outlier_label.reordered, outlier_label),
        outlier_label2 = coalesce(outlier_label2.reordered, outlier_label2)
    ) %>%
    select(-ref.reordered, -outlier_label.reordered, -outlier_label2.reordered)%>%
  mutate(count = sapply(strsplit(outliers, ";"), length))

# Filter outliers by label
outliers_control <- ehdn_results_annotated %>%
    filter(outlier_label %in% c("control"))

outliers_cases <- ehdn_results_annotated %>%
    filter(outlier_label == "case")

outliers_mixed <- ehdn_results_annotated %>%
    filter(outlier_label == "mixed")

# Annotate data with additional information
filter_and_annotate <- function(df) {
    df <- df %>% 
        mutate(
            motif_length = nchar(motif),
            expansion_length = end - start,
            Expansion_Double_Motif = expansion_length >= 2 * motif_length,
            type = ifelse(nchar(motif) <= 6, 'STR', 'VNTR')
        ) 
    return(df)
}

outliers_1_case <- ehdn_results_annotated %>% filter(count == 1 & outlier_label == "case")

outliers_cases <- filter_and_annotate(outliers_cases)
outliers_control <- filter_and_annotate(outliers_control)
outliers_mixed <- filter_and_annotate(outliers_mixed)
outliers_1_case <- filter_and_annotate(outliers_1_case)

# Function to plot motif length comparison
plot_motif_length_comparison <- function(dfs, cohort_names, colors, Title, filename) {
    calculate_motif_lengths <- function(df) {
        motif_lengths <- nchar(as.character(df$motif))
        motif_length_counts <- table(motif_lengths)
        motif_length_counts_df <- data.frame(length = as.numeric(names(motif_length_counts)), count = as.numeric(motif_length_counts))
        return(motif_length_counts_df)
    }

    motif_length_counts_list <- list()
    for (i in seq_along(dfs)) {
        motif_length_counts_df <- calculate_motif_lengths(dfs[[i]])
        motif_length_counts_df <- mutate(motif_length_counts_df, cohort = cohort_names[i])
        motif_length_counts_list[[i]] <- motif_length_counts_df
    }

    ehdn_motif_length_counts_combined <- bind_rows(motif_length_counts_list)
    motif_length_comparison_plot <- ggplot(ehdn_motif_length_counts_combined, aes(x = length, y = count, fill = cohort)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = Title, x = "Motif Length", y = "Count") +
        scale_fill_manual(values = colors) +
        theme_minimal() +
        theme(text = element_text(family = "Times", color = "black"))
    
    ggsave(filename, plot = motif_length_comparison_plot)
    print(motif_length_comparison_plot)
}

# Function to plot motif length comparison as percentages
plot_motif_length_comparison_percentage <- function(dfs, cohort_names, colors, Title, filename) {
    calculate_percentage <- function(df) { 
        motif_lengths_3 <- nchar(as.character(df$motif))
        motif_length_counts_3 <- table(motif_lengths_3)
        motif_length_counts_df_3 <- data.frame(length = as.numeric(names(motif_length_counts_3)), count = as.numeric(motif_length_counts_3))
        motif_length_counts_df_3$percentage <- (motif_length_counts_df_3$count / sum(motif_length_counts_df_3$count)) * 100
        return(motif_length_counts_df_3)
    }

    percentage_dfs <- list()
    for (i in seq_along(dfs)) {
        df_percentage <- calculate_percentage(dfs[[i]])
        df_percentage <- mutate(df_percentage, cohort = cohort_names[i])
        percentage_dfs[[i]] <- df_percentage
    }

    ehdn_motiflength_percentage_combined <- bind_rows(percentage_dfs)
    motif_length_comparison_plot_percentage <- ggplot(ehdn_motiflength_percentage_combined, aes(x = length, y = percentage, fill = cohort)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = Title, x = "Nucleotide Length", y = "Percentage") +
        scale_fill_manual(values = colors) +
        theme_minimal() +
        theme(text = element_text(family = "Times", color = "black"))
    
    ggsave(filename, plot = motif_length_comparison_plot_percentage)
    print(motif_length_comparison_plot_percentage)
}

# Function to plot nucleotide comparison
plot_nucleotide_comparison <- function(dfs, cohort_names, colors, Title, filename) {
    calculate_nucleotide <- function(df) {
        nucleotide <- table(unlist(strsplit(as.character(df$motif), "")))
        nucleotide_counts <- data.frame(base = names(nucleotide), count = as.numeric(nucleotide))
        nucleotide_counts_df <- nucleotide_counts[!nucleotide_counts$base %in% "N", ]
        return(nucleotide_counts_df)
    }

    nucleotide_counts_list <- list()
    for (i in seq_along(dfs)) {
        nucleotide_counts_df <- calculate_nucleotide(dfs[[i]])
        nucleotide_counts_df <- mutate(nucleotide_counts_df, cohort = cohort_names[i])
        nucleotide_counts_list[[i]] <- nucleotide_counts_df
    }

    nucleotide_counts_combined <- bind_rows(nucleotide_counts_list)
    nucleotide_comparison_plot <- ggplot(nucleotide_counts_combined, aes(x = base, y = count, fill = cohort)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = Title, x = "Nucleotide", y = "Count") +
        scale_fill_manual(values = colors) +
        theme_minimal() +
        theme(text = element_text(family = "Times", color = "black"))
    
    ggsave(filename, plot = nucleotide_comparison_plot)
    print(nucleotide_comparison_plot)
}

# Function to plot nucleotide comparison as percentages
plot_nucleotide_comparison_percentage <- function(dfs, cohort_names, colors, Title, filename) {
    calculate_nucleotide <- function(df) {
        nucleotide <- table(unlist(strsplit(as.character(df$motif), "")))
        nucleotide_counts <- data.frame(base = names(nucleotide), count = as.numeric(nucleotide))
        nucleotide_counts_df <- nucleotide_counts[!nucleotide_counts$base %in% "N", ]
        nucleotide_counts_df$percentage <- (nucleotide_counts_df$count / sum(nucleotide_counts_df$count)) * 100
        return(nucleotide_counts_df)
    }

    nucleotide_counts_list <- list()
    for (i in seq_along(dfs)) {
        nucleotide_counts_df <- calculate_nucleotide(dfs[[i]])
        nucleotide_counts_df <- mutate(nucleotide_counts_df, cohort = cohort_names[i])
        nucleotide_counts_list[[i]] <- nucleotide_counts_df
    }

    nucleotide_counts_combined <- bind_rows(nucleotide_counts_list)
    nucleotide_comparison_plot_percentage <- ggplot(nucleotide_counts_combined, aes(x = base, y = percentage, fill = cohort)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = Title, x = "Nucleotide", y = "Percentage") +
        scale_fill_manual(values = colors) +
        theme_minimal() +
        theme(text = element_text(family = "Times", color = "black"))
    
    ggsave(filename, plot = nucleotide_comparison_plot_percentage)
    print(nucleotide_comparison_plot_percentage)
}

# Function to plot region comparison
plot_region_comparison <- function(dfs, cohort_names, colors, Title, filename) {
    for (i in seq_along(dfs)) {
        dfs[[i]] <- mutate(dfs[[i]], cohort = cohort_names[i])
    }

    region_combined <- bind_rows(dfs)
    region_comparison_plot <- ggplot(region_combined, aes(x = region, fill = cohort)) +
        geom_bar(stat = "count", position = "dodge") +
        labs(title = Title, x = "Region", y = "Count") +
        scale_fill_manual(values = colors) +
        theme_minimal(base_size = 5) +
        theme(text = element_text(family = "Times", color = "black"))
    
    ggsave(filename, plot = region_comparison_plot)
    print(region_comparison_plot)
}

# Function to plot region comparison as percentages
plot_region_comparison_percentage <- function(dfs, cohort_names, colors, Title, filename) {
    calculate_region <- function(df) {
        count <- table(df$region)
        count_1 <- data.frame(region = names(count), count = as.vector(count))
        count_1$percentage <- (count_1$count / sum(count_1$count)) * 100
        return(count_1)
    }

    region_list <- list()
    for (i in seq_along(dfs)) {
        region_df <- calculate_region(dfs[[i]])
        region_df <- mutate(region_df, cohort = cohort_names[i])
        region_list[[i]] <- region_df
    }

    region_combined <- bind_rows(region_list)
    region_comparison_plot <- ggplot(region_combined, aes(x = region, y = percentage, fill = cohort)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = Title, x = "Region", y = "Percentage of Tandem Repeats") +
        scale_fill_manual(values = colors) +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(text = element_text(family = "Times", color = "black"))
    
    ggsave(filename, plot = region_comparison_plot)
    print(region_comparison_plot)
}

# Function to plot prominent motif comparison
plot_motif_comparison <- function(dfs, cohort_names, colors, Title, filename) {
    calculate_motif <- function(df) {
        motif_counts <- table(df$motif)
        motif <- data.frame(motif = names(motif_counts), count = as.vector(motif_counts))
        sorted_motif <- motif[order(-motif$count), ]
        top_10_motifs <- sorted_motif[1:10, ]
        return(top_10_motifs)
    }

    motif_list <- list()
    for (i in seq_along(dfs)) {
        top_10_motifs <- calculate_motif(dfs[[i]])
        top_10_motifs <- mutate(top_10_motifs, cohort = cohort_names[i])
        motif_list[[i]] <- top_10_motifs
    }

    motif_combined <- bind_rows(motif_list)
    motif_comparison_plot <- ggplot(motif_combined, aes(x = motif, y = count, fill = cohort)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = Title, x = "Motif", y = "Count") +
        scale_fill_manual(values = colors) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(text = element_text(family = "Times", color = "black"))
    
    ggsave(filename, plot = motif_comparison_plot)
    print(motif_comparison_plot)
}

# Function to plot prominent motif comparison as percentages
plot_motif_comparison_percentage <- function(dfs, cohort_names, colors, Title, filename) {
    calculate_motif <- function(df) {
        motif_counts <- table(df$motif)
        motif <- data.frame(motif = names(motif_counts), count = as.vector(motif_counts))
        sorted_motif <- motif[order(-motif$count), ]
        total_count <- sum(motif$count)
        sorted_motif$percentage <- (sorted_motif$count / total_count) * 100    
        top_10_motifs <- sorted_motif[1:10, ]
        return(top_10_motifs)
    }

    motif_list <- list()
    for (i in seq_along(dfs)) {
        top_10_motifs <- calculate_motif(dfs[[i]])
        top_10_motifs <- mutate(top_10_motifs, cohort = cohort_names[i])
        motif_list[[i]] <- top_10_motifs
    }

    motif_combined <- bind_rows(motif_list)
    motif_comparison_percentage_plot <- ggplot(motif_combined, aes(x = motif, y = percentage, fill = cohort)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = Title, x = "Motif", y = "Percentage") +
        scale_fill_manual(values = colors) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(text = element_text(family = "Times", color = "black"))
    
    ggsave(filename, plot = motif_comparison_percentage_plot)
    print(motif_comparison_percentage_plot)
}

# Function to plot chromosome comparison
plot_cromosome_comparison <- function(dfs, cohort_names, colors, Title, filename) {
    calculate_chromosome <- function(df) {
        cromosome <- table(unlist(as.character(df$contig), ""))
        cromosome <- data.frame(base = names(cromosome), count = as.numeric(cromosome))
        return(cromosome)
    }

    cromosome_counts_list <- list()
    for (i in seq_along(dfs)) {
        cromosome <- calculate_chromosome(dfs[[i]])
        cromosome <- mutate(cromosome, cohort = cohort_names[i])
        cromosome_counts_list[[i]] <- cromosome
    }

    cromosome_counts_combined <- bind_rows(cromosome_counts_list)
    cromosome_comparison_plot <- ggplot(cromosome_counts_combined, aes(x = base, y = count, fill = cohort)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = Title, x = "Motif Length", y = "Count") +
        scale_fill_manual(values = colors) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(text = element_text(family = "Times", color = "black"))
    
    ggsave(filename, plot = cromosome_comparison_plot)
    print(cromosome_comparison_plot)
}

# Function to plot chromosome comparison as percentages
plot_cromosome_comparison_percentage <- function(dfs, cohort_names, colors, Title, filename) {
    calculate_chromosome <- function(df) {
        cromosome <- table(unlist(as.character(df$contig), ""))
        cromosome <- data.frame(base = names(cromosome), count = as.numeric(cromosome))
        total_count <- sum(cromosome$count)
        cromosome$percentage <- (cromosome$count / total_count) * 100
        return(cromosome)
    }

    cromosome_counts_list <- list()
    for (i in seq_along(dfs)) {
        cromosome <- calculate_chromosome(dfs[[i]])
        cromosome <- mutate(cromosome, cohort = cohort_names[i])
        cromosome_counts_list[[i]] <- cromosome
    }

    cromosome_counts_combined <- bind_rows(cromosome_counts_list)
    cromosome_comparison_percentage_plot <- ggplot(cromosome_counts_combined, aes(x = base, y = percentage, fill = cohort)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = Title, x = "Motif Length", y = "Count") +
        scale_fill_manual(values = colors) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(text = element_text(family = "Times", color = "black"))
    
    ggsave(filename, plot = cromosome_comparison_percentage_plot)
    print(cromosome_comparison_percentage_plot)
}

# Function to plot nucleotide comparison for STRs and VNTRs
plot_nucleotide_comparison_str_vntr <- function(dfs, cohort_names, colors, Title, filename) {
    calculate_nucleotide <- function(df) {
        nucleotide <- table(unlist(strsplit(as.character(df$motif), "")))
        nucleotide_counts <- data.frame(base = names(nucleotide), count = as.numeric(nucleotide))
        nucleotide_counts_df <- nucleotide_counts[!nucleotide_counts$base %in% "N", ]
        return(nucleotide_counts_df)
    }

    nucleotide_counts_list <- list()
    for (i in seq_along(dfs)) {
        nucleotide_counts_df <- calculate_nucleotide(dfs[[i]])
        nucleotide_counts_df <- mutate(nucleotide_counts_df, cohort = cohort_names[i])
        nucleotide_counts_list[[i]] <- nucleotide_counts_df
    }

    nucleotide_counts_combined <- bind_rows(nucleotide_counts_list)
    nucleotide_comparison_plot <- ggplot(nucleotide_counts_combined, aes(x = base, y = count, fill = cohort)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = Title, x = "Nucleotide", y = "Count") +
        scale_fill_manual(values = colors) +
        theme_minimal() +
        theme(text = element_text(family = "Times", color = "black"))
    
    ggsave(filename, plot = nucleotide_comparison_plot)
    print(nucleotide_comparison_plot)
}

# Function to plot nucleotide comparison for STRs and VNTRs as percentages
plot_nucleotide_comparison_str_vntr_percentage <- function(dfs, cohort_names, colors, Title, filename) {
    calculate_nucleotide <- function(df) {
        nucleotide <- table(unlist(strsplit(as.character(df$motif), "")))
        nucleotide_counts <- data.frame(base = names(nucleotide), count = as.numeric(nucleotide))
        nucleotide_counts_df <- nucleotide_counts[!nucleotide_counts$base %in% "N", ]
        nucleotide_counts_df$percentage <- (nucleotide_counts_df$count / sum(nucleotide_counts_df$count)) * 100
        return(nucleotide_counts_df)
    }

    nucleotide_counts_list <- list()
    for (i in seq_along(dfs)) {
        nucleotide_counts_df <- calculate_nucleotide(dfs[[i]])
        nucleotide_counts_df <- mutate(nucleotide_counts_df, cohort = cohort_names[i])
        nucleotide_counts_list[[i]] <- nucleotide_counts_df
    }

    nucleotide_counts_combined <- bind_rows(nucleotide_counts_list)
    nucleotide_comparison_percentage_plot <- ggplot(nucleotide_counts_combined, aes(x = base, y = percentage, fill = cohort)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = Title, x = "Nucleotide", y = "Percentage") +
        scale_fill_manual(values = colors) +
        theme_minimal() +
        theme(text = element_text(family = "Times", color = "black"))
    
    ggsave(filename, plot = nucleotide_comparison_percentage_plot)
    print(nucleotide_comparison_percentage_plot)
}

# Function to plot expansion length comparison with scatter plot
plot_expansion_length_comparison_scatter <- function(dfs, cohort_names, colors, Title, filename) {
    expansion_length_list <- list()
    for (i in seq_along(dfs)) {
        dfs[[i]] <- mutate(dfs[[i]], cohort = cohort_names[i])
        expansion_length_list[[i]] <- dfs[[i]]
    }

    expansion_length_combined <- bind_rows(expansion_length_list)
    expansion_length_comparison_scatter_plot <- ggplot(expansion_length_combined, aes(x = motif_length, y = expansion_length, color = cohort)) +
        geom_point(alpha = 0.5) +
        labs(title = Title, x = "Repeat Motif Length", y = "Repeat Expansion Length") +
        scale_color_manual(values = colors) +
        theme_minimal() +
        theme(text = element_text(family = "Times", color = "black"))
    
    ggsave(filename, plot = expansion_length_comparison_scatter_plot)
    print(expansion_length_comparison_scatter_plot)
}

# Function to plot expansion length comparison with box plot
plot_expansion_length_comparison_box <- function(dfs, cohort_names, colors, Title, filename) {
    expansion_length_list <- list()
    for (i in seq_along(dfs)) {
        dfs[[i]] <- mutate(dfs[[i]], cohort = cohort_names[i])
        expansion_length_list[[i]] <- dfs[[i]]
    }

    expansion_length_combined <- bind_rows(expansion_length_list)
    expansion_length_comparison_box_plot <- ggplot(expansion_length_combined, aes(x = factor(motif_length), y = expansion_length, fill = cohort)) +
        geom_boxplot() +
        labs(title = Title, x = "Repeat Motif Length", y = "Repeat Expansion Length") +
        scale_fill_manual(values = colors) +
        theme_minimal() +
        theme(text = element_text(family = "Times", color = "black"))
    
    ggsave(filename, plot = expansion_length_comparison_box_plot)
    print(expansion_length_comparison_box_plot)
}

# Apply the functions to the specified dataframes and plot the results
plot_motif_length_comparison(
    list(outliers_cases, outliers_control, outliers_mixed, outliers_1_case),
    c("case", "control", "mixed", "1_case"),
    c("purple", "blue", "green", "magenta"),
    "Motif Length Comparison of Outliers",
    "motif_length_comparison.png"
)

plot_motif_length_comparison_percentage(
    list(outliers_cases, outliers_control, outliers_mixed, outliers_1_case),
    c("case", "control", "mixed", "1_case"),
    c("purple", "blue", "green", "magenta"),
    "Motif Length Comparison of Outliers (Percentage)",
    "motif_length_comparison_percentage.png"
)

plot_nucleotide_comparison(
    list(outliers_cases, outliers_control, outliers_mixed, outliers_1_case),
    c("case", "control", "mixed", "1_case"),
    c("purple", "blue", "green", "magenta"),
    "Nucleotide Comparison of Outliers",
    "nucleotide_comparison.png"
)

plot_nucleotide_comparison_percentage(
    list(outliers_cases, outliers_control, outliers_mixed, outliers_1_case),
    c("case", "control", "mixed", "1_case"),
    c("purple", "blue", "green", "magenta"),
    "Nucleotide Comparison of Outliers (Percentage)",
    "nucleotide_comparison_percentage.png"
)

plot_region_comparison(
    list(outliers_cases, outliers_control, outliers_mixed, outliers_1_case),
    c("case", "control", "mixed", "1_case"),
    c("purple", "blue", "green", "magenta"),
    "Region Comparison of Outliers",
    "region_comparison.png"
)

plot_region_comparison_percentage(
    list(outliers_cases, outliers_control, outliers_mixed, outliers_1_case),
    c("case", "control", "mixed", "1_case"),
    c("purple", "blue", "green", "magenta"),
    "Region Comparison of Outliers (Percentage)",
    "region_comparison_percentage.png"
)

plot_motif_comparison(
    list(outliers_cases, outliers_control, outliers_mixed, outliers_1_case),
    c("case", "control", "mixed", "1_case"),
    c("purple", "blue", "green", "magenta"),
    "Motif Comparison of Outliers",
    "motif_comparison.png"
)

plot_motif_comparison_percentage(
    list(outliers_cases, outliers_control, outliers_mixed, outliers_1_case),
    c("case", "control", "mixed", "1_case"),
    c("purple", "blue", "green", "magenta"),
    "Motif Comparison of Outliers (Percentage)",
    "motif_comparison_percentage.png"
)

plot_cromosome_comparison(
    list(outliers_cases, outliers_control, outliers_mixed, outliers_1_case),
    c("case", "control", "mixed", "1_case"),
    c("purple", "blue", "green", "magenta"),
    "Chromosome Comparison of Outliers",
    "chromosome_comparison.png"
)

plot_cromosome_comparison_percentage(
    list(outliers_cases, outliers_control, outliers_mixed, outliers_1_case),
    c("case", "control", "mixed", "1_case"),
    c("purple", "blue", "green", "magenta"),
    "Chromosome Comparison of Outliers (percentage)",
    "chromosome_comparison_percentage.png"
)

plot_nucleotide_comparison_str_vntr(
    list(outliers_cases %>% filter(type == "STR"), outliers_control %>% filter(type == "STR"), outliers_mixed %>% filter(type == "STR"), outliers_1_case %>% filter(type == "STR")),
    c("case", "control", "mixed", "1_case"),
    c("purple", "blue", "green", "magenta"),
    "Nucleotide Comparison of STRs",
    "nucleotide_comparison_str.png"
)

plot_nucleotide_comparison_str_vntr(
    list(outliers_cases %>% filter(type == "VNTR"), outliers_control %>% filter(type == "VNTR"), outliers_mixed %>% filter(type == "VNTR"), outliers_1_case %>% filter(type == "VNTR")),
    c("case", "control", "mixed", "1_case"),
    c("purple", "blue", "green", "magenta"),
    "Nucleotide Comparison of VNTRs",
    "nucleotide_comparison_vntr.png"
)

plot_nucleotide_comparison_str_vntr_percentage(
    list(outliers_cases %>% filter(type == "STR"), outliers_control %>% filter(type == "STR"), outliers_mixed %>% filter(type == "STR"), outliers_1_case %>% filter(type == "STR")),
    c("case", "control", "mixed", "1_case"),
    c("purple", "blue", "green", "magenta"),
    "Nucleotide Comparison of STRs (Percentage)",
    "nucleotide_comparison_str_percentage.png"
)

plot_nucleotide_comparison_str_vntr_percentage(
    list(outliers_cases %>% filter(type == "VNTR"), outliers_control %>% filter(type == "VNTR"), outliers_mixed %>% filter(type == "VNTR"), outliers_1_case %>% filter(type == "VNTR")),
    c("case", "control", "mixed", "1_case"),
    c("purple", "blue", "green", "magenta"),
    "Nucleotide Comparison of VNTRs (Percentage)",
    "nucleotide_comparison_vntr_percentage.png"
)

plot_expansion_length_comparison_scatter(
    list(outliers_cases, outliers_control, outliers_mixed, outliers_1_case),
    c("case", "control", "mixed", "1_case"),
    c("purple", "blue", "green", "magenta"),
    "Expansion Length Comparison of Outliers (Scatter Plot)",
    "expansion_length_comparison_scatter.png"
)

plot_expansion_length_comparison_box(
    list(outliers_cases, outliers_control, outliers_mixed, outliers_1_case),
    c("case", "control", "mixed", "1_case"),
    c("purple", "blue", "green", "magenta"),
    "Expansion Length Comparison of Outliers (Box Plot)",
    "expansion_length_comparison_box.png"
)



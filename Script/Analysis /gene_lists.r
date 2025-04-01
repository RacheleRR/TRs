#BRAIN FILTER 
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(kableExtra)
library(tidyverse)
library(readxl)
library(broom)

# load 
brain_gene_consensus_filtered_consensus_no_pitular <- read.delim("~/brain_gene_list/consensus/brain_gene_consensus_filtered_consensus_no_pitular_no_unique.tsv")
brain_gene_consensus_ntm_consensus_no_pitular <- read.delim("~/brain_gene_list/consensus/brain_gene_consensus_ntm_consensus_no_pitular_no_unique.tsv")
ehdn_results_annotated <- read.delim("~/EHDN_DBSCAN_correct/Result/results_dbscan_with_QC_NOUHR/AFTER_DBSCAN/ehdn_DBSCAN_annotated.tsv")
ehdn_results_reorder <- read.delim("~/EHDN_DBSCAN_correct/Result/results_dbscan_with_QC_NOUHR/AFTER_DBSCAN/ehdn_DBSCAN_reorder.tsv")

#get labels 
ehdn_results_annotated <- ehdn_results_annotated %>%
  left_join(ehdn_results_reorder %>% select(repeatID, ref, outlier_label, outlier_label2), by = "repeatID", suffix = c("", ".reordered")) %>%
  mutate(
    ref = coalesce(ref.reordered, ref),
    outlier_label = coalesce(outlier_label.reordered, outlier_label),
    outlier_label2 = coalesce(outlier_label2.reordered, outlier_label2)
  ) %>%
  select(-ref.reordered, -outlier_label.reordered, -outlier_label2.reordered) %>%
  mutate(count = sapply(strsplit(outliers, ";"), length))

#prepare data for filtering split genes 
ehdn_results_annotated_split <- ehdn_results_annotated %>%
  rename("Gene" = gene) %>%
  mutate(Gene = gsub('\\([^\\)]+\\)', '', Gene)) %>%
  separate_rows(Gene, sep = ',') %>%
  mutate(CpG = ifelse(grepl("CG|GC|CGC|CGG|CCG|GCC|GGC|GCG", motif), 1, 0))

#keep only genic regions
ehdn_genic <- ehdn_results_annotated_split %>% filter(region != "intergenic")
ehdn_exonic <- ehdn_results_annotated_split %>% filter(region == "exonic")
ehdn_genic_pure <- ehdn_results_annotated_split %>% filter(region != "intergenic" & region != "upstream" & region != "downstream" & region != "ncRNA_intronic" & region != "ncRNA_exonic" & region != "ncRNA_splicing" & region != "ncRNA_UTR3" & region != "ncRNA_UTR5" & region != "upstream;downstream" & region != "downstream;upstream")

#keep only in Case and count 1 and exonic or genic 
ehdn_exonic_Case1 <- ehdn_exonic %>% filter(outlier_label == "Case", count == 1)
ehdn_genic_Case1 <- ehdn_genic %>% filter(outlier_label == "Case", count == 1)
ehdn_genic_Case <- ehdn_genic %>% filter(outlier_label == "Case")
ehdn_exonic_Case <- ehdn_exonic %>% filter(outlier_label == "Case")
ehdn_genic_Control <- ehdn_genic %>% filter(outlier_label == "Control")
ehdn_exonic_Control <- ehdn_exonic %>% filter(outlier_label == "Control")
ehdn_genic_mixed <- ehdn_genic %>% filter(outlier_label == "mixed")
ehdn_exonic_mixed <- ehdn_exonic %>% filter(outlier_label == "mixed")
ehdn_genic_pure_Case1 <- ehdn_genic_pure %>% filter(outlier_label == "Case", count == 1)
ehdn_genic_pure_Case <- ehdn_genic_pure %>% filter(outlier_label == "Case")
ehdn_genic_pure_Control <- ehdn_genic_pure %>% filter(outlier_label == "Control")
ehdn_genic_pure_mixed <- ehdn_genic_pure %>% filter(outlier_label == "mixed")

# filter by Case and Case count 1 and CPG 1 (doesn't matter the region)
ehdn_Case1_cpg1 <- ehdn_results_annotated_split %>% filter(outlier_label == "Case", count == 1, CpG == 1)
ehdn_Case_cpg <- ehdn_results_annotated_split %>% filter(outlier_label == "Case", CpG == 1)

#BRAIN LISTS
ehdn <- ehdn_results_annotated_split
#separate it into Case, Control, mixed 
genes_brain_ntmp <- unique(brain_gene_consensus_ntm_consensus_no_pitular$Gene.name)
genes_brain_filter <- unique(brain_gene_consensus_filtered_consensus_no_pitular$Gene.name)
common_brain_ntmp <- intersect(ehdn$Gene, genes_brain_ntmp)
common_brain_filter <- intersect(ehdn$Gene, genes_brain_filter)

# Aggregate Tissue information
brain_gene_consensus_ntm_consensus_no_pitular <- brain_gene_consensus_ntm_consensus_no_pitular %>%
  group_by(Gene.name) %>%
  summarise(Tissue = paste(unique(Tissue), collapse = ","))
brain_gene_consensus_filtered_consensus_no_pitular <- brain_gene_consensus_filtered_consensus_no_pitular %>%
  group_by(Gene.name) %>%
  summarise(Tissue = paste(unique(Tissue), collapse = ","))

#filter gene list by Case, Control, mixed and Case1
Case_brain_ntmp <- ehdn %>% filter(Gene %in% common_brain_ntmp, outlier_label == "Case") %>%
  left_join(brain_gene_consensus_ntm_consensus_no_pitular, by = c("Gene" = "Gene.name"))
Control_brain_ntmp <- ehdn %>% filter(Gene %in% common_brain_ntmp, outlier_label == "Control") %>%
  left_join(brain_gene_consensus_ntm_consensus_no_pitular, by = c("Gene" = "Gene.name"))
mixed_brain_ntmp <- ehdn %>% filter(Gene %in% common_brain_ntmp, outlier_label == "mixed") %>%
  left_join(brain_gene_consensus_ntm_consensus_no_pitular, by = c("Gene" = "Gene.name"))
Case_brain_filter <- ehdn %>% filter(Gene %in% common_brain_filter, outlier_label == "Case") %>%
  left_join(brain_gene_consensus_filtered_consensus_no_pitular, by = c("Gene" = "Gene.name"))
Control_brain_filter <- ehdn %>% filter(Gene %in% common_brain_filter, outlier_label == "Control") %>%
  left_join(brain_gene_consensus_filtered_consensus_no_pitular, by = c("Gene" = "Gene.name"))
mixed_brain_filter <- ehdn %>% filter(Gene %in% common_brain_filter, outlier_label == "mixed") %>%
  left_join(brain_gene_consensus_filtered_consensus_no_pitular, by = c("Gene" = "Gene.name"))
Case_brain_ntmp_Case1 <- ehdn %>% filter(Gene %in% common_brain_ntmp, outlier_label == "Case", count == 1) %>%
  left_join(brain_gene_consensus_ntm_consensus_no_pitular, by = c("Gene" = "Gene.name"))
Case_brain_filter_Case1 <- ehdn %>% filter(Gene %in% common_brain_filter, outlier_label == "Case", count == 1) %>%
  left_join(brain_gene_consensus_filtered_consensus_no_pitular, by = c("Gene" = "Gene.name"))


#filter gene list by Case and exonic and genic 
Case_brain_ntmp_exonic <- Case_brain_ntmp %>% filter(region == "exonic")
Case_brain_ntmp_genic <- Case_brain_ntmp %>% filter(region != "intergenic")
Case_brain_filter_exonic <- Case_brain_filter %>% filter(region == "exonic")
Case_brain_filter_genic <- Case_brain_filter %>% filter(region != "intergenic")
Case_brain_ntmp_Case1_exonic <- Case_brain_ntmp_Case1 %>% filter(region == "exonic")
Case_brain_ntmp_Case1_genic <- Case_brain_ntmp_Case1 %>% filter(region != "intergenic")
Case_brain_filter_Case1_exonic <- Case_brain_filter_Case1 %>% filter(region == "exonic")
Case_brain_filter_Case1_genic <- Case_brain_filter_Case1 %>% filter(region != "intergenic")
Case_brain_filter_genic_pure <- Case_brain_filter %>% filter(region != "intergenic" & region != "upstream" & region != "downstream" & region != "ncRNA_intronic" & region != "ncRNA_exonic" & region != "ncRNA_splicing" & region != "ncRNA_UTR3" & region != "ncRNA_UTR5" & region != "upstream;downstream" & region != "downstream;upstream")
Case_brain_ntmp_genic_pure <- Case_brain_ntmp %>% filter(region != "intergenic" & region != "upstream" & region != "downstream" & region != "ncRNA_intronic" & region != "ncRNA_exonic" & region != "ncRNA_splicing" & region != "ncRNA_UTR3" & region != "ncRNA_UTR5" & region != "upstream;downstream" & region != "downstream;upstream")
Case_brain_ntmp_Case1_genic_pure <- Case_brain_ntmp_Case1 %>% filter(region != "intergenic" & region != "upstream" & region != "downstream" & region != "ncRNA_intronic" & region != "ncRNA_exonic" & region != "ncRNA_splicing" & region != "ncRNA_UTR3" & region != "ncRNA_UTR5" & region != "upstream;downstream" & region != "downstream;upstream")
Case_brain_filter_Case1_genic_pure <- Case_brain_filter_Case1 %>% filter(region != "intergenic" & region != "upstream" & region != "downstream" & region != "ncRNA_intronic" & region != "ncRNA_exonic" & region != "ncRNA_splicing" & region != "ncRNA_UTR3" & region != "ncRNA_UTR5" & region != "upstream;downstream" & region != "downstream;upstream")

#SAVE
setwd("/home/rachele/EHDN_DBSCAN_correct/Result/results_dbscan_with_QC_NOUHR/AFTER_DBSCAN/gene_lists")

# Create a summary dataframe with the count of unique genes in each dataframe
library(kableExtra)
library(dplyr)

# Create the summary dataframe
summary_df <- data.frame(
  DataFrame = c(
    "ehdn_results_annotated_split", "ehdn_results_annotated", "ehdn_genic", "ehdn_exonic", "ehdn_genic_pure",
    "ehdn_exonic_Case1", "ehdn_genic_Case1", "ehdn_genic_Case", "ehdn_exonic_Case", "ehdn_genic_Control",
    "ehdn_exonic_Control", "ehdn_genic_mixed", "ehdn_exonic_mixed", "ehdn_genic_pure_Case", "ehdn_genic_pure_Control",
    "ehdn_genic_pure_mixed", "ehdn_genic_pure_Case1", "ehdn_Case1_cpg1", "ehdn_Case_cpg", "Case_brain_ntmp",
    "Control_brain_ntmp", "mixed_brain_ntmp", "Case_brain_filter", "Control_brain_filter", "mixed_brain_filter",
    "Case_brain_ntmp_Case1", "Case_brain_filter_Case1", "Case_brain_ntmp_exonic", "Case_brain_ntmp_genic",
    "Case_brain_filter_exonic", "Case_brain_filter_genic", "Case_brain_ntmp_Case1_exonic", "Case_brain_ntmp_Case1_genic",
    "Case_brain_filter_Case1_exonic", "Case_brain_filter_Case1_genic", "Case_brain_filter_genic_pure", "Case_brain_ntmp_genic_pure",
    "Case_brain_ntmp_Case1_genic_pure", "Case_brain_filter_Case1_genic_pure"
  ),
  UniqueGenes = c(
    length(unique(ehdn_results_annotated_split$Gene)), length(unique(ehdn_results_annotated$Gene)), length(unique(ehdn_genic$Gene)),
    length(unique(ehdn_exonic$Gene)), length(unique(ehdn_genic_pure$Gene)), length(unique(ehdn_exonic_Case1$Gene)),
    length(unique(ehdn_genic_Case1$Gene)), length(unique(ehdn_genic_Case$Gene)), length(unique(ehdn_exonic_Case$Gene)),
    length(unique(ehdn_genic_Control$Gene)), length(unique(ehdn_exonic_Control$Gene)), length(unique(ehdn_genic_mixed$Gene)),
    length(unique(ehdn_exonic_mixed$Gene)), length(unique(ehdn_genic_pure_Case$Gene)), length(unique(ehdn_genic_pure_Control$Gene)),
    length(unique(ehdn_genic_pure_mixed$Gene)), length(unique(ehdn_genic_pure_Case1$Gene)), length(unique(ehdn_Case1_cpg1$Gene)),
    length(unique(ehdn_Case_cpg$Gene)), length(unique(Case_brain_ntmp$Gene)), length(unique(Control_brain_ntmp$Gene)),
    length(unique(mixed_brain_ntmp$Gene)), length(unique(Case_brain_filter$Gene)), length(unique(Control_brain_filter$Gene)),
    length(unique(mixed_brain_filter$Gene)), length(unique(Case_brain_ntmp_Case1$Gene)), length(unique(Case_brain_filter_Case1$Gene)),
    length(unique(Case_brain_ntmp_exonic$Gene)), length(unique(Case_brain_ntmp_genic$Gene)), length(unique(Case_brain_filter_exonic$Gene)),
    length(unique(Case_brain_filter_genic$Gene)), length(unique(Case_brain_ntmp_Case1_exonic$Gene)), length(unique(Case_brain_ntmp_Case1_genic$Gene)),
    length(unique(Case_brain_filter_Case1_exonic$Gene)), length(unique(Case_brain_filter_Case1_genic$Gene)), length(unique(Case_brain_filter_genic_pure$Gene)),
    length(unique(Case_brain_ntmp_genic_pure$Gene)), length(unique(Case_brain_ntmp_Case1_genic_pure$Gene)), length(unique(Case_brain_filter_Case1_genic_pure$Gene))
  )
)

# Use kableExtra to create a styled table
summary_table <- summary_df %>%
  kable("html", caption = "Summary of Unique Genes Across DataFrames") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
  scroll_box(width = "100%", height = "500px")

# Save the summary table as an HTML file
save_kable(summary_table, file = "summary_unique_genes.html")

# Save the summary dataframe as a tsv file
write.table(summary_df, file = "summary_unique_genes.tsv", sep = "\t", row.names = FALSE)

# Save the summary dataframe as a tsv file
write.table(summary_df, file = "summary_unique_genes.tsv", sep = "\t", row.names = FALSE)

summary_df_1 <- data.frame(
  DataFrame = c(
    "ehdn_results_annotated_split", "ehdn_results_annotated", "ehdn_genic", "ehdn_exonic", "ehdn_genic_pure",
    "ehdn_exonic_Case1", "ehdn_genic_Case1", "ehdn_genic_Case", "ehdn_exonic_Case", "ehdn_genic_Control",
    "ehdn_exonic_Control", "ehdn_genic_mixed", "ehdn_exonic_mixed", "ehdn_genic_pure_Case", "ehdn_genic_pure_Control",
    "ehdn_genic_pure_mixed", "ehdn_genic_pure_Case1", "ehdn_Case1_cpg1", "ehdn_Case_cpg"
  ),
  UniqueGenes = c(
    length(unique(ehdn_results_annotated_split$Gene)), length(unique(ehdn_results_annotated$Gene)), length(unique(ehdn_genic$Gene)),
    length(unique(ehdn_exonic$Gene)), length(unique(ehdn_genic_pure$Gene)), length(unique(ehdn_exonic_Case1$Gene)),
    length(unique(ehdn_genic_Case1$Gene)), length(unique(ehdn_genic_Case$Gene)), length(unique(ehdn_exonic_Case$Gene)),
    length(unique(ehdn_genic_Control$Gene)), length(unique(ehdn_exonic_Control$Gene)), length(unique(ehdn_genic_mixed$Gene)),
    length(unique(ehdn_exonic_mixed$Gene)), length(unique(ehdn_genic_pure_Case$Gene)), length(unique(ehdn_genic_pure_Control$Gene)),
    length(unique(ehdn_genic_pure_mixed$Gene)), length(unique(ehdn_genic_pure_Case1$Gene)), length(unique(ehdn_Case1_cpg1$Gene)),
    length(unique(ehdn_Case_cpg$Gene))
  )
)

summary_df_2 <- data.frame(
  DataFrame = c(
    "Case_brain_ntmp",
    "Control_brain_ntmp", "mixed_brain_ntmp", "Case_brain_filter", "Control_brain_filter", "mixed_brain_filter",
    "Case_brain_ntmp_Case1", "Case_brain_filter_Case1"
  ),
  UniqueGenes = c(
    length(unique(Case_brain_ntmp$Gene)), length(unique(Control_brain_ntmp$Gene)),
    length(unique(mixed_brain_ntmp$Gene)), length(unique(Case_brain_filter$Gene)), length(unique(Control_brain_filter$Gene)),
    length(unique(mixed_brain_filter$Gene)), length(unique(Case_brain_ntmp_Case1$Gene)), length(unique(Case_brain_filter_Case1$Gene))  
  )
)

summary_df_3 <- data.frame(
  DataFrame = c(
    "Case_brain_ntmp_exonic", "Case_brain_ntmp_genic",
    "Case_brain_filter_exonic", "Case_brain_filter_genic", "Case_brain_ntmp_Case1_exonic", "Case_brain_ntmp_Case1_genic",
    "Case_brain_filter_Case1_exonic", "Case_brain_filter_Case1_genic", "Case_brain_filter_genic_pure", "Case_brain_ntmp_genic_pure",
    "Case_brain_ntmp_Case1_genic_pure", "Case_brain_filter_Case1_genic_pure"
  ),
  UniqueGenes = c(
    length(unique(Case_brain_ntmp_exonic$Gene)), length(unique(Case_brain_ntmp_genic$Gene)), length(unique(Case_brain_filter_exonic$Gene)),
    length(unique(Case_brain_filter_genic$Gene)), length(unique(Case_brain_ntmp_Case1_exonic$Gene)), length(unique(Case_brain_ntmp_Case1_genic$Gene)),
    length(unique(Case_brain_filter_Case1_exonic$Gene)), length(unique(Case_brain_filter_Case1_genic$Gene)), length(unique(Case_brain_filter_genic_pure$Gene)),
    length(unique(Case_brain_ntmp_genic_pure$Gene)), length(unique(Case_brain_ntmp_Case1_genic_pure$Gene)), length(unique(Case_brain_filter_Case1_genic_pure$Gene))
  )
)

#save them all as tsv files, sep \t and row name false , dont give the path
write.table(ehdn_results_annotated_split, file = "ehdn_results_annotated_split.tsv", sep = "\t", row.names = FALSE)
write.table(ehdn_results_annotated, file = "ehdn_results_annotated_with_labels.tsv", sep = "\t", row.names = FALSE)

write.table(ehdn_genic, file = "ehdn_genic.tsv", sep = "\t", row.names = FALSE)
write.table(ehdn_exonic, file = "ehdn_exonic.tsv", sep = "\t", row.names = FALSE)
write.table(ehdn_genic_pure, file = "ehdn_genic_pure.tsv", sep = "\t", row.names = FALSE)

write.table(ehdn_exonic_Case1, file = "ehdn_exonic_Case1.tsv", sep = "\t", row.names = FALSE)
write.table(ehdn_genic_Case1, file = "ehdn_genic_Case1.tsv", sep = "\t", row.names = FALSE)
write.table(ehdn_genic_Case, file = "ehdn_genic_Case.tsv", sep = "\t", row.names = FALSE)
write.table(ehdn_exonic_Case, file = "ehdn_exonic_Case.tsv", sep = "\t", row.names = FALSE)
write.table(ehdn_genic_Control, file = "ehdn_genic_Control.tsv", sep = "\t", row.names = FALSE)     
write.table(ehdn_exonic_Control, file = "ehdn_exonic_Control.tsv", sep = "\t", row.names = FALSE)
write.table(ehdn_genic_mixed, file = "ehdn_genic_mixed.tsv", sep = "\t", row.names = FALSE)
write.table(ehdn_exonic_mixed, file = "ehdn_exonic_mixed.tsv", sep = "\t", row.names = FALSE)
write.table(ehdn_genic_pure_Case, file = "ehdn_genic_pure_Case.tsv", sep = "\t", row.names = FALSE)
write.table(ehdn_genic_pure_Control, file = "ehdn_genic_pure_Control.tsv", sep = "\t", row.names = FALSE)
write.table(ehdn_genic_pure_mixed, file = "ehdn_genic_pure_mixed.tsv", sep = "\t", row.names = FALSE)
write.table(ehdn_genic_pure_Case1, file = "ehdn_genic_pure_Case1.tsv", sep = "\t", row.names = FALSE)


write.table(ehdn_Case1_cpg1, file = "ehdn_Case1_cpg1.tsv", sep = "\t", row.names = FALSE)
write.table(ehdn_Case_cpg, file = "ehdn_Case_cpg.tsv", sep = "\t", row.names = FALSE)

write.table(Case_brain_ntmp, file = "Case_brain_ntmp.tsv", sep = "\t", row.names = FALSE)
write.table(Control_brain_ntmp, file = "Control_brain_ntmp.tsv", sep = "\t", row.names = FALSE)
write.table(mixed_brain_ntmp, file = "mixed_brain_ntmp.tsv", sep = "\t", row.names = FALSE)
write.table(Case_brain_filter, file = "Case_brain_filter.tsv", sep = "\t", row.names = FALSE)
write.table(Control_brain_filter, file = "Control_brain_filter.tsv", sep = "\t", row.names = FALSE)
write.table(mixed_brain_filter, file = "mixed_brain_filter.tsv", sep = "\t", row.names = FALSE)
write.table(Case_brain_ntmp_Case1, file = "Case_brain_ntmp_Case1.tsv", sep = "\t", row.names = FALSE)
write.table(Case_brain_filter_Case1, file = "Case_brain_filter_Case1.tsv", sep = "\t", row.names = FALSE)


write.table(Case_brain_ntmp_exonic, file = "Case_brain_ntmp_exonic.tsv", sep = "\t", row.names = FALSE)
write.table(Case_brain_ntmp_genic, file = "Case_brain_ntmp_genic.tsv", sep = "\t", row.names = FALSE)
write.table(Case_brain_filter_exonic, file = "Case_brain_filter_exonic.tsv", sep = "\t", row.names = FALSE)
write.table(Case_brain_filter_genic, file = "Case_brain_filter_genic.tsv", sep = "\t", row.names = FALSE)
write.table(Case_brain_ntmp_Case1_exonic, file = "Case_brain_ntmp_Case1_exonic.tsv", sep = "\t", row.names = FALSE)
write.table(Case_brain_ntmp_Case1_genic, file = "Case_brain_ntmp_Case1_genic.tsv", sep = "\t", row.names = FALSE)
write.table(Case_brain_filter_Case1_exonic, file = "Case_brain_filter_Case1_exonic.tsv", sep = "\t", row.names = FALSE)
write.table(Case_brain_filter_Case1_genic, file = "Case_brain_filter_Case1_genic.tsv", sep = "\t", row.names = FALSE)
write.table(Case_brain_filter_genic_pure, file = "Case_brain_filter_genic_pure.tsv", sep = "\t", row.names = FALSE)
write.table(Case_brain_ntmp_genic_pure, file = "Case_brain_ntmp_genic_pure.tsv", sep = "\t", row.names = FALSE)
write.table(Case_brain_ntmp_Case1_genic_pure, file = "Case_brain_ntmp_Case1_genic_pure.tsv", sep = "\t", row.names = FALSE)
write.table(Case_brain_filter_Case1_genic_pure, file = "Case_brain_filter_Case1_genic_pure.tsv", sep = "\t", row.names = FALSE)

# Function to create and save Excel files
save_to_excel <- function(data_list, summary_sheet, file_name) {
  final_list <- c(list(Summary = summary_sheet), data_list)
  write.xlsx(final_list, file = file_name, rowNames = FALSE)
}

# Create a list of data frames
data_list_1 <- list(
  summary_numbers = summary_df_1,  
  ehdn_genic = ehdn_genic,
  ehdn_exonic = ehdn_exonic,
  ehdn_genic_pure = ehdn_genic_pure,
  ehdn_exonic_Case1 = ehdn_exonic_Case1,
  ehdn_genic_Case1 = ehdn_genic_Case1,
  ehdn_genic_Case = ehdn_genic_Case,
  ehdn_exonic_Case = ehdn_exonic_Case,
  ehdn_genic_Control = ehdn_genic_Control,
  ehdn_exonic_Control = ehdn_exonic_Control,
  ehdn_genic_mixed = ehdn_genic_mixed,
  ehdn_exonic_mixed = ehdn_exonic_mixed,
  ehdn_genic_pure_Case = ehdn_genic_pure_Case,
  ehdn_genic_pure_Control = ehdn_genic_pure_Control,
  ehdn_genic_pure_mixed = ehdn_genic_pure_mixed,
  ehdn_genic_pure_Case1 = ehdn_genic_pure_Case1,
  ehdn_Case1_cpg1 = ehdn_Case1_cpg1,
  ehdn_Case_cpg = ehdn_Case_cpg
)

# Create a summary data frame for the first sheet
summary_sheet_1 <- data.frame(
  DataFrame_Name = names(data_list_1),
  Description = c(
    "Summary of unique genes in each dataframe.",
    "Genic regions from EHDN analysis.",
    "Exonic regions from EHDN analysis.",
    "Pure genic regions from EHDN analysis.",
    "Exonic regions from EHDN analysis (Case1 subset).",
    "Genic regions from EHDN analysis (Case1 subset).",
    "Genic regions from EHDN analysis (Case subset).",
    "Exonic regions from EHDN analysis (Case subset).",
    "Genic regions from EHDN analysis (Control subset).",
    "Exonic regions from EHDN analysis (Control subset).",
    "Genic regions from EHDN analysis (mixed subset).",
    "Exonic regions from EHDN analysis (mixed subset).",
    "Pure genic regions from EHDN analysis (Case subset).",
    "Pure genic regions from EHDN analysis (Control subset).",
    "Pure genic regions from EHDN analysis (mixed subset).",
    "Pure genic regions from EHDN analysis (Case1 subset).",
    "CpG-related data from EHDN analysis (Case1 subset).",
    "CpG-related data from EHDN analysis (Case subset)."
  )
)

# Save the first set of data frames to an Excel file
save_to_excel(data_list_1, summary_sheet_1, "ehdn_data.xlsx")

# Create a list of data frames
data_list_2 <- list(
  summary_numbers = summary_df_2,
  Case_brain_ntmp = Case_brain_ntmp,
  Control_brain_ntmp = Control_brain_ntmp,
  mixed_brain_ntmp = mixed_brain_ntmp,
  Case_brain_filter = Case_brain_filter,
  Control_brain_filter = Control_brain_filter,
  mixed_brain_filter = mixed_brain_filter,
  Case_brain_ntmp_Case1 = Case_brain_ntmp_Case1,
  Case_brain_filter_Case1 = Case_brain_filter_Case1
)

# Create a summary data frame for the first sheet
summary_sheet_2 <- data.frame(
  DataFrame_Name = names(data_list_2),
  Description = c(
    "Summary of unique genes in each brain dataframe.",
    "Case brain data (no filtering).",
    "Control brain data (no filtering).",
    "Mixed brain data (no filtering).",
    "Case brain data (filtered).",
    "Control brain data (filtered).",
    "Mixed brain data (filtered).",
    "Case brain data (Case1 subset, no filtering).",
    "Case brain data (Case1 subset, filtered)."
  )
)

# Save the second set of data frames to an Excel file
save_to_excel(data_list_2, summary_sheet_2, "brain_data.xlsx")


# Create a list of data frames
data_list_3 <- list(
  summary_numbers = summary_df_3,  
  Case_brain_ntmp_exonic = Case_brain_ntmp_exonic,
  Case_brain_ntmp_genic = Case_brain_ntmp_genic,
  Case_brain_filter_exonic = Case_brain_filter_exonic,
  Case_brain_filter_genic = Case_brain_filter_genic,
  Case_brain_ntmp_Case1_exonic = Case_brain_ntmp_Case1_exonic,
  Case_brain_ntmp_Case1_genic = Case_brain_ntmp_Case1_genic,
  Case_brain_filter_Case1_exonic = Case_brain_filter_Case1_exonic,
  Case_brain_filter_Case1_genic = Case_brain_filter_Case1_genic,
  Case_brain_filter_genic_pure = Case_brain_filter_genic_pure,
  Case_brain_ntmp_genic_pure = Case_brain_ntmp_genic_pure,
  Case1_brain_ntmp_genic_pure = Case_brain_ntmp_Case1_genic_pure,
  Case1_brain_filt_genic_pure = Case_brain_filter_Case1_genic_pure
)

# Create a summary data frame for the first sheet
summary_sheet_3 <- data.frame(
  DataFrame_Name = names(data_list_3),
  Description = c(
    "Summary of unique genes in each brain and region dataframe.",
    "Exonic regions from Case brain (no filtering).",
    "Genic regions from Case brain (no filtering).",
    "Exonic regions from Case brain (filtered).",
    "Genic regions from Case brain (filtered).",
    "Exonic regions from Case brain (Case1 subset, no filtering).",
    "Genic regions from Case brain (Case1 subset, no filtering).",
    "Exonic regions from Case brain (Case1 subset, filtered).",
    "Genic regions from Case brain (Case1 subset, filtered).",
    "Pure genic regions from Case brain (filtered).",
    "Pure genic regions from Case brain (no filtering).",
    "Pure genic regions from Case brain (Case1 subset, no filtering).",
    "Pure genic regions from Case brain (Case1 subset, filtered)."
  )
)

# Save the third set of data frames to an Excel file
save_to_excel(data_list_3, summary_sheet_3, "Case_brain_data_region.xlsx")








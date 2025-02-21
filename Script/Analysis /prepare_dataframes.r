#Process EHDN DBSCAN

# prepare the ehdn annotated
library(dplyr)
library(tidyr)
library(tidyverse)

# Read the reordered and annotated TSV files
ehdn_results_reordered <- read.delim("/home/rachele/ehdn_DBSCAN_reorder.tsv", stringsAsFactors = FALSE)
ehdn_results_annotated <- read.delim("/home/rachele/ehdn_DBSCAN_annotated.tsv", stringsAsFactors = FALSE)


ehdn_results_annotated <- ehdn_results_annotated %>%
  left_join(ehdn_results_reordered %>% select(repeatID, ref, outlier_label, outlier_label2), by = "repeatID", suffix = c("", ".reordered")) %>%
  mutate(
    ref = coalesce(ref.reordered, ref),
    outlier_label = coalesce(outlier_label.reordered, outlier_label),
    outlier_label2 = coalesce(outlier_label2.reordered, outlier_label2)
  ) %>%
  select(-ref.reordered, -outlier_label.reordered, -outlier_label2.reordered)


#analysis 1 outlier in 1 individual ( ORA) general


# Function to filter and annotate the data frame
filter_and_annotate <- function(df) {
    df <- df %>%
        rename("Gene" = gene) %>%
        mutate(Gene = gsub('\\([^\\)]+\\)', '', Gene)) %>%
        separate_rows(Gene, sep = ',') %>%
        mutate(
            motif_length = nchar(motif),
            expansion_length = end - start,
            Expansion_Double_Motif = expansion_length >= 2 * motif_length,
            T_count = str_count(motif, 'T'),
            C_count = str_count(motif, 'C'),
            G_count = str_count(motif, 'G'),
            A_count = str_count(motif, 'A'),
            type = ifelse(nchar(motif) <= 6, 'STR', 'VNTR')
        ) %>%
        mutate(
            CG_percentage = ((C_count + G_count) / motif_length) * 100,
            AT_percentage = ((A_count + T_count) / motif_length) * 100
        ) %>%
        mutate(
            T_percentage = (T_count / motif_length) * 100,
            C_percentage = (C_count / motif_length) * 100,
            G_percentage = (G_count / motif_length) * 100,
            A_percentage = (A_count / motif_length) * 100
        )

    return(df)
}

ehdn_modified <- filter_and_annotate(ehdn_results_annotated)

ehdn_modified <- ehdn_modified %>%
  mutate(count = sapply(strsplit(outliers, ";"), length))

outliers_1_individual_case <- ehdn_modified %>%
  filter(count == 1 & outlier_label == "case")

# analysis 1 outlier in 1 indiviual exonic

Outlier_1_exonic_case <- ehdn_modified %>% filter(region == "exonic" & count == 1 & outlier_label =="case")

# TSS + junction
outliers_control_mixed <- ehdn_results_annotated %>%
    filter(outlier_label %in% c("control", "mixed"))


outliers_cases <- ehdn_results_annotated %>%
    filter(count > 1  & outlier_label == "case")




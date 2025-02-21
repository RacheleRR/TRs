library(tidyr)
library(dplyr)

manifest_noUHRNA <- read.table("~/TRs/manifest_noUHRNA.tsv", sep = "\t", header = TRUE)
Ehdn_dbscan_noUHR_noQC <- read.table("~/EHDN_DBSCAN_correct/Result/results_dbscan_without_QC_NOUHR/EHdn.expansions.2025-02-18.tsv", sep = "\t" , header = T)


# Define the input file with the outlier values
input_file <- "Ehdn_dbscan_noUHR_noQC"
# Define the manifest file
 manifest_file <- "manifest_noUHRNA"
# 
# # Load the manifest dataframe
 manifest_df <- get(manifest_file)

# Load the input dataframe
input_df <- get(input_file)

#check if they are case or control 



#WAY 1 ==> just tell me if pur control, case or mized

# Function to check the label for the outlier values
check_outlier_label <- function(outlier_value, manifest_df) {
    # Split the outlier values by ';' and trim whitespace
    outlier_values <- strsplit(outlier_value, ";")[[1]]
    
    # Check if any of the outlier values exist in the manifest dataframe
    labels <- sapply(outlier_values, function(val) {
        # Return the label (case or control) from the manifest if the value matches
        if (val %in% manifest_df$V1) {
            return(manifest_df$V2[manifest_df$V1 == val])
        } else {
            return(NA) # Return NA if no match
        }
    })
    
    # If there's more than one label, return the unique ones
    unique_labels <- unique(labels)
    
    # If there's both "case" and "control", return "mixed"; otherwise return the label
    if (length(unique_labels) > 1) {
        return("mixed")
    } else if (length(unique_labels) == 1) {
        return(unique_labels)
    } else {
        return(NA) # If no label found
    }
}

# Apply the function to the 'outliers' column and create a new column 'outlier_label'
input_df$outlier_label <- sapply(input_df$outliers, check_outlier_label, manifest_df = manifest_df)


#WAY2 ==> tell me exactly what label 

# Function to get the labels (case or control) for the outlier values
get_outlier_labels <- function(outlier_value, manifest_df) {
    # Split the outlier values by ';' into a vector
    outlier_values <- strsplit(outlier_value, ";")[[1]]
    
    # Retrieve corresponding labels from the manifest dataframe
    labels <- sapply(outlier_values, function(val) {
        # Find the label in the manifest dataframe (either "case" or "control")
        if (val %in% manifest_df$V1) {
            return(manifest_df$V2[manifest_df$V1 == val])
        } else {
            return(NA) # Return NA if no match is found
        }
    })
    
    # Concatenate the labels into a single string, separated by commas
    return(paste(labels, collapse = ", "))
}

# Apply the function to the 'outliers' column and create a new column 'outlier_label2'
input_df$outlier_label2 <- sapply(input_df$outliers, get_outlier_labels, manifest_df = manifest_df)

# Assign the modified dataframe back to the original name
assign(input_file, input_df)

# Count the occurrences of 'case', 'control', and 'mixed'
case_count <- sum(grepl("case", input_df$outlier_label) & !grepl("control", input_df$outlier_label))
control_count <- sum(grepl("control", input_df$outlier_label) & !grepl("case", input_df$outlier_label))
mixed_count <- sum(grepl("mixed", input_df$outlier_label) & !grepl("control", input_df$outlier_label))

# Display the counts
cat("Number of 'case' labels:", case_count, "\n")
cat("Number of 'control' labels:", control_count, "\n")
cat("Number of 'mixed' labels:", mixed_count, "\n")


#clean it up to use for annovar 
ehdn_results_reordered <- input_df %>%
rename(contig = chr) %>%
select(contig, start, end, everything())


# Write the reordered data to a new TSV file
write.table(ehdn_results_reordered, "/home/rachele/EHDN_DBSCAN_correct/Result/results_dbscan_without_QC_NOUHR/AFTER_DBSCAN/ehdn_DBSCAN_reorder.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
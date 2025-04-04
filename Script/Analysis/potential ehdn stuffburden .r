# Analysis of the EHDN data


#FROM HERE
annotation.genic <- ehdn_results_annotated$repeatID[
    !is.na(ehdn_results_annotated$gene) & 
        ehdn_results_annotated$region %in% c("exonic", "exonic;splicing", "upstream", "downstream", "UTR3", "UTR5", "splicing", "intronic","ncRNA_intronic","ncRNA_exonic","ncRNA_splicing")
]

annotation.intergenic <- ehdn_results_annotated$repeatID[
    ehdn_results_annotated$region == "intergenic"
]


# Step 1: Define the subset of loci that are:
#   - On chromosomes 1-22 (`contig %in% paste0("chr", 1:22)`)
#   - In genic regions (`repeatID %in% annotation.genic`)
subset.genic <- ehdn_results_annotated$contig %in% paste0("chr", 1:22) & 
                ehdn_results_annotated$repeatID %in% annotation.genic

# Step 2: Extract the `outliers` (sample lists) for these loci
ehdn.genic.outliers <- ehdn_results_annotated$outliers[subset.genic]

# Step 3: Collapse all outliers into a single string and split into individual samples
ehdn.genic.samples <- strsplit(paste(ehdn.genic.outliers, collapse = ";"), ";")[[1]]

# Step 4: Count the number of expansions per sample
ehdn.genic.counts <- data.frame(table(ehdn.genic.samples), stringsAsFactors = FALSE)
names(ehdn.genic.counts) <- c("Sample", "EHdnGenotypeGenic")


# define subset for intergenic 

subset.intergenic <- ehdn_results_annotated$contig %in% paste0("chr", 1:22) & 
                ehdn_results_annotated$repeatID %in% annotation.intergenic
ehdn.intergenic.outliers <- ehdn_results_annotated$outliers[subset.intergenic ]
ehdn.intergenic.samples <- strsplit(paste(ehdn.intergenic.outliers, collapse = ";"), ";")[[1]]
ehdn.intergenic.counts <- data.frame(table(ehdn.intergenic.samples), stringsAsFactors = FALSE)
names(ehdn.intergenic.counts) <- c("Sample", "EHdnGenotypeIntergenic")


merge.data <- merge(ehdn.intergenic.counts, ehdn.genic.counts, by = "Sample", all = T)
names(merge.data) <- c("Sample", "intergenic", "genic")
merge.data[is.na(merge.data)] <- 0


#get status from the manifest file
manifest_df <- read.csv("~/manifest_SCZdetail.csv")

#if the same sample id get the status from the manifest file
names(manifest_df) <- c("Status", "Sex", "Age", "Sample")
merge.data <- merge(merge.data, manifest_df, by = "Sample", all = T)

#if its NA in the genic or intergenic  then put 0 
merge.data[is.na(merge.data)] <- 0


lm <- glm(Status ~ genic ,data=merge.data , family = binomial)
summary(lm)
exp(lm$coefficients["genic"])
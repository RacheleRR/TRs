# candidate GENES

#? LOAD 
 #load library 
        library(readr)
        library(data.table)
        # library(tidyverse)
        library(dplyr)
        library(gprofiler2)
        library(httr)
        library(readxl)
        library(tidyr)
        library(stringr)
        setwd("/home/rachele/EHDN_DBSCAN_correct/Result/results_dbscan_with_QC_NOUHR/AFTER_DBSCAN/")
 # load my data   
        ehdn  <- read.delim("outliers_1_case_split.tsv")
#        ehdn  <- read.delim("outliers_1_case_exonic_split.tsv")

ehdn <- ehdn %>%   mutate(CpG = ifelse(grepl("CG|GC|CGC|CGG|CCG|GCC|GGC|GCG", motif), 1, 0))


 #LOAD GENESETS
        SCHEMA <- read_csv("/home/rachele/Documents/old/geneset_confrontation/SCHEMA.csv")
        GWAS_120 <- read.delim("~/GWAS_120.csv")       
      
#? clean GENESETS
 # CONVERT GENE NAMES 

    # SCHEMA
        convert_col <- gconvert(query = SCHEMA$Gene, organism = "hsapiens", target = "ENSG", mthreshold = Inf, filter_na = TRUE)
        SCHEMA <- merge(SCHEMA, convert_col[, c("target", "name")], by.x = "Gene", by.y = "target", all.x = TRUE)
        SCHEMA <- SCHEMA %>% select(Gene, name, everything())


 # General clean 

                                  SCHEMA <- SCHEMA%>%subset(select = -Gene)%>% rename('Gene' = name )%>% mutate(Gene=gsub(' ','',Gene)) %>%setNames(make.names(colnames(.)))  
                                  # #SCHEMAp.V meno 0.05 
                                  SCHEMA_pVal <- SCHEMA[SCHEMA$P.meta <= 0.05, ]     
                                  GWAS_120 <- GWAS_120 %>% rename('Gene' = GENE_name )
                                  
#? UNIFY ALL THE DATAFRAMES 

        # Extract "Gene" column from each dataframe
        genes_schema_pval <- SCHEMA_pVal$Gene
        genes_gwas <- GWAS_120$Gene


        # Ensure uniqueness in the "Gene" column for each list
        genes_schema_pval <- unique(genes_schema_pval)
        genes_gwas <- unique(genes_gwas)

#? intersect with ehdn
        common_ehdn_schizophrenia <- intersect(ehdn$Gene, genes_schema_pval)
        common_ehdn_gwas <- intersect(ehdn$Gene, genes_gwas)


common_ehdn_schizophrenia <-ehdn  %>% filter(Gene %in% genes_schema_pval)
common_ehdn_gwas <- ehdn  %>% filter(Gene %in% genes_gwas)

write.csv(common_ehdn_schizophrenia, file = "common_ehdn_schizophrenia.csv", row.names = FALSE)
write.csv(common_ehdn_gwas, file = "common_ehdn_gwas.csv", row.names = FALSE)

# filter for no intergenic and yes CPG 1
common_ehdn_schizophrenia_cpg <- common_ehdn_schizophrenia %>% filter(CpG == 1 & region != "intergenic")
# filter for no intergenic and yes CPG 1
common_ehdn_gwas_cpg <- common_ehdn_gwas %>% filter(CpG == 1 & region != "intergenic")

write.csv(common_ehdn_schizophrenia_cpg, file = "common_ehdn_schizophrenia_cpg.csv", row.names = FALSE)
write.csv(common_ehdn_gwas_cpg, file = "common_ehdn_gwas_cpg.csv", row.names = FALSE)



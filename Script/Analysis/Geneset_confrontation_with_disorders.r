#! 13) Geneset confrontation with disorders 
#TODO: To asses overlap between my Geneset and Pathogenic Genesets ( no statistics)
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
        setwd("/home/rachele/EHDN_DBSCAN_correct/Result/results_dbscan_without_QC_NOUHR/AFTER_DBSCAN/")
 # load my data   
        ehdn  <- read.delim("outliers_1_case_split.tsv")
#        ehdn  <- read.delim("outliers_1_case_exonic_split.tsv")



 #LOAD GENESETS
        SCHEMA <- read_csv("/home/rachele/Documents/old/geneset_confrontation/SCHEMA.csv")
        GWAS_120 <- read.delim("~/GWAS_120.csv")       
        BipEx_Bipolar <- read_csv("/home/rachele/Documents/old/geneset_confrontation/Bipolar_Disorder.csv")
        SFARI <- read_csv("/home/rachele/Documents/old/geneset_confrontation/SFARI.csv")

#? clean GENESETS
 # CONVERT GENE NAMES 
    # BipEx_Bipolar
        convert_col <- gconvert(query = BipEx_Bipolar$Gene, organism = "hsapiens", target = "ENSG", mthreshold = Inf, filter_na = TRUE)
        BipEx_Bipolar <- merge(BipEx_Bipolar, convert_col[, c("target", "name")], by.x = "Gene", by.y = "target", all.x = TRUE)
        BipEx_Bipolar <- BipEx_Bipolar %>% select(Gene, name, everything())

    # SCHEMA
        convert_col <- gconvert(query = SCHEMA$Gene, organism = "hsapiens", target = "ENSG", mthreshold = Inf, filter_na = TRUE)
        SCHEMA <- merge(SCHEMA, convert_col[, c("target", "name")], by.x = "Gene", by.y = "target", all.x = TRUE)
        SCHEMA <- SCHEMA %>% select(Gene, name, everything())


 # General clean 
                                  BipEx_Bipolar <- BipEx_Bipolar%>%subset(select = -Gene)%>% rename('Gene' = name )%>% mutate(Gene=gsub(' ','',Gene))%>%  setNames(make.names(colnames(.))) 

                                BipEx_Bipolar_p_val_PTV <- BipEx_Bipolar[!is.na(BipEx_Bipolar$PTV.Fisher.p.val) & BipEx_Bipolar$PTV.Fisher.p.val <= 0.05,]
                                BipEx_Bipolar_p_val_Missense <- BipEx_Bipolar[!is.na(BipEx_Bipolar$Damaging.Missense.Fisher.p.val) & BipEx_Bipolar$Damaging.Missense.Fisher.p.val <= 0.05,]


                                  SCHEMA <- SCHEMA%>%subset(select = -Gene)%>% rename('Gene' = name )%>% mutate(Gene=gsub(' ','',Gene)) %>%setNames(make.names(colnames(.)))  
                                  # # SCHEMA OR +1
                                  SCHEMA_OR <- SCHEMA[SCHEMA$OR..Class.I. >= 1, ]
                                  # #SCHEMAp.V meno 0.05 
                                  SCHEMA_pVal <- SCHEMA[SCHEMA$P.meta <= 0.05, ]
                                  SFARI <- SFARI %>%subset(select = -c(status, `ensembl-id`))%>%rename("Gene"=`gene-symbol`, "gene_name"=`gene-name`)%>% mutate(Gene=gsub(' ','',Gene))
                                  SFARI_syndromic <- SFARI 
                                  SFARI_non_syndromic_lower_3 <- SFARI[SFARI$`gene-score` < 3, ]
                                  
                                  GWAS_120 <- GWAS_120 %>% rename('Gene' = GENE_name )
                                  
#? UNIFY ALL THE DATAFRAMES 

        # Extract "Gene" column from each dataframe
        genes_schema_pval <- SCHEMA_pVal$Gene
        genes_bipolar <- c(BipEx_Bipolar_p_val_PTV$Gene, BipEx_Bipolar_p_val_Missense$Gene)
        genes_sfari_non_syndromic <- SFARI_non_syndromic_lower_3$Gene
        genes_gwas <- GWAS_120$Gene
        genes_schema_or <- SCHEMA_OR$Gene
        genes_sfari_syndromic <- SFARI_syndromic$Gene

        # Ensure uniqueness in the "Gene" column for each list
        genes_schema_pval <- unique(genes_schema_pval)
        genes_bipolar <- unique(genes_bipolar)
        genes_sfari_non_syndromic <- unique(genes_sfari_non_syndromic)
        genes_gwas <- unique(genes_gwas)
        genes_schema_or <- unique(genes_schema_or)
        genes_sfari_syndromic <- unique(genes_sfari_syndromic)

#? intersect with ehdn
        common_ehdn_schizophrenia <- intersect(ehdn$Gene, genes_schema_pval)
        common_ehdn_bipolar <- intersect(ehdn$Gene, genes_bipolar)
        common_ehdn_autism <- intersect(ehdn$Gene, genes_sfari_non_syndromic)
        common_ehdn_gwas <- intersect(ehdn$Gene, genes_gwas)
        common_ehdn_schema_or <- intersect(ehdn$Gene, genes_schema_or)
        common_ehdn_sfari_syndromic <- intersect(ehdn$Gene, genes_sfari_syndromic)

#? VENNDIAGRAM 
        library(VennDiagram)
        library(RColorBrewer)
        
        # Define custom colors
        venn_colors <- c("#FF6347", "#4682B4", "#32CD32", "#FFD700", "#8A2BE2")  # Example colors: Tomato, SteelBlue, LimeGreen, Gold, BlueViolet
        
        # First Venn Diagram
        x1 <- list(SCHEMA_pVal = common_ehdn_schizophrenia, Bipolar = common_ehdn_bipolar, SFARI_non_syndromic = common_ehdn_autism, GWAS = common_ehdn_gwas, SCHEMA_OR = common_ehdn_schema_or)
        
        # Save the first Venn diagram as a PNG file
        png("Venn_Schema_Bipolar_SFARI_GWAS_OR.png", width = 800, height = 800)
        grid.newpage()
        venn.plot1 <- venn.diagram(
          x1,
          filename = NULL,
          col = venn_colors,
          fill = venn_colors,
          alpha = 0.5,
          cex = 1,
          cat.cex = 1,
          cat.col = venn_colors,
          main = "Venn Diagram 1: SCHEMA_pVal, Bipolar, SFARI_non_syndromic, GWAS, SCHEMA_OR"
        )
        grid.draw(venn.plot1)
        dev.off()  # Close the graphics device
        
        # Second Venn Diagram
        common_ehdn_gwas_schema_pval <- unique(c(common_ehdn_gwas, common_ehdn_schizophrenia))
        x2 <- list(GWAS_SCHEMA_pVal = common_ehdn_gwas_schema_pval, Bipolar = common_ehdn_bipolar, SFARI_syndromic = common_ehdn_sfari_syndromic, SFARI_non_syndromic = common_ehdn_autism, SCHEMA_OR = common_ehdn_schema_or)
        
        # Save the second Venn diagram as a PNG file
        png("Venn_GWAS_SFARI_SCHEMABipolar.png", width = 800, height = 800)
        grid.newpage()
        venn.plot2 <- venn.diagram(
          x2,
          filename = NULL,
          col = venn_colors,
          fill = venn_colors,
          alpha = 0.5,
          cex = 1,
          cat.cex = 1,
          cat.col = venn_colors,
          main = "Venn Diagram 2: GWAS_SCHEMA_pVal, Bipolar, SFARI_syndromic, SFARI_non_syndromic, SCHEMA_OR"
        )
        grid.draw(venn.plot2)
        dev.off()  # Close the graphics device
        
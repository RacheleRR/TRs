# TRs
It illustrates the process of identifying outliers using EHDN and DBSCAN, followed by other analysis 

# FIND OUTLIERS USING EHDN AND DBSCAN
====================================

## Overview
This project demonstrates the process of identifying outliers using EHDN (Expansion Hunter Denovo) and DBSCAN (Density-Based Spatial Clustering of Applications with Noise). The steps outlined below detail how the data is processed and analyzed.


## Objectives
- Identify outliers in the dataset.
- Utilize EHDN for profiling and merging calls.
- Apply DBSCAN for clustering and outlier detection.



## Prerequisites
- [EHDN](https://github.com/Illumina/ExpansionHunterDenovo)
- R
- Python
- Fazzar2020Script()
- BTlib ()

## Steps

1. **Apply Profiling Command**
    - Run the `ehdn_new_calls.sh` script to perform the profiling step using EHDN.

2. **Perform Quality Control**
    - Execute the QC.sh script to perform quality control checks on the data and identifes and removes samples as outliers from anifesto file 

3. **Combine Samples&Regions and perform DBSCAN**
   -Run the Combine_RUN_DBSCAN.sh to aggregate sample data, merge and compare regions and apply DBSCAN to cluster and detect outliers. 

4. **Annotate**
   -Run Annotate_C_C_prep_ANNOVAR.r to prep for the application of annotation script using Annovar and apply annotate_annovar.sh

5. **Clean Dataframe and prepare it for future analysis**
   -Run prepare_dataframes.r

6. **Correlation analysis**
  - To investigate the association between genomic features and TRE. Run Correlation_in_domains.r

8. **Proximity Analysis**
   - To assess spatial relationships with regulatory element and TRE. Run TSS_junction.r
     
10. **Statistical Analysis**
    - To perform fisher, kruskal and wilcoxon analysis. RUn #
   



## Commands run 

**Run profiling**
./ehdn_new_calls.sh

**Perform QC** 
 ./QC.sh -i /home/rachele/EHdn/new_calls/ -m /home/rachele/manifest_noUHRNA.tsv  -o /home/rachele/EHDN_DBSCAN_correct/Result/QC

**Combine regions and samples and apply DBSCAN** 

./Combine_RUN_DBSCAN.sh -m /home/rachele/EHDN_DBSCAN_correct/Result/QC/filtered_manifest.tsv

./Combine_RUN_DBSCAN.sh -m /home/rachele/manifest_noUHRNA.tsv 




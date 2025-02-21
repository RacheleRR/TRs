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
- DBSCAN



## Steps

1. **Apply Profiling Command**
    - Run the `ehdn_new_calls.sh` script to perform the profiling step using EHDN.

2. **Perform Quality Control**
    - Execute the QC.sh script to perform quality control checks on the data and identifes and removes samples as outliers from anifesto file 

3. **Combine Samples +Regions and DBSCAN**
   -Run the Combine_RUN_DBSCAN.sh to aggregate sample data, merge and compare regions and apply DBSCAN to cluster and detect outliers. 

4. **



## Commands ruN 

# Run profiling
./ehdn_new_calls.sh

# Perform QC 
 ./QC.sh -i /home/rachele/EHdn/new_calls/ -m /home/rachele/manifest_noUHRNA.tsv  -o /home/rachele/EHDN_DBSCAN_correct/Result/QC


# Combine regions and samples and apply DBSCAN 

./Combine_RUN_DBSCAN.sh -m /home/rachele/EHDN_DBSCAN_correct/Result/QC/filtered_manifest.tsv

./Combine_RUN_DBSCAN.sh -m /home/rachele/manifest_noUHRNA.tsv 




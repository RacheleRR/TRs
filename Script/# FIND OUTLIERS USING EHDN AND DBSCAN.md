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
    - Execute the `count_calls.R` script to perform quality control checks on the data.

3. **Remove Outliers**
    - Identify and remove samples flagged as outliers from the manifesto file.

4. **Combine Samples**
    - Run the `combine_counts.py` script to aggregate sample data. (Script provided by Fazal.)

5. **Combine Regions**
    - Use the `compare_anchored_irrs.py` script to merge and compare regions. (Script provided by Fazal.)

6. **Find Outliers Using DBSCAN**
    - Apply DBSCAN for clustering and outlier detection by running the `dbscan_only.R` script.

7. **See if Case,Control or Mixed**
    - Apply script to either see overview or have detailed desciption by running 'Case_orControl.R' script. 



## Example Commands

```bash
# Run profiling
./ehdn_new_calls.sh

# Perform quality control
Rscript count_calls.R

# Combine samples+Combine regions+ DBSCAN
./Combine_RUN_DBSCAN.sh 


```






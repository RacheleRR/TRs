# 🧬 TRs: Finding Outliers with EHDN & DBSCAN

Welcome to **TRs**! This project guides you through the process of detecting outliers in genomic data using **EHDN (Expansion Hunter Denovo)** and **DBSCAN (Density-Based Spatial Clustering of Applications with Noise)**. 

---

## 🌟 Overview
This repository demonstrates how to:
- Identify outliers in genomic datasets
- Use **EHDN** for profiling and merging calls
- Apply **DBSCAN** for clustering and outlier detection
- Perform downstream statistical and proximity analyses

---

## 🎯 Objectives
✔️ Detect outliers in genomic data  
✔️ Process and refine data using **EHDN**  
✔️ Cluster and analyze data using **DBSCAN**  
✔️ Annotate, clean, and correlate results for further analysis  

---

## 🛠️ Prerequisites
Make sure you have the following installed:
- [EHDN](https://github.com/Illumina/ExpansionHunterDenovo) 🧬
- **R** 📊
- **Python** 🐍
- **Fazzar2020Script()**
- **BTlib ()**

---

## 🔍 Steps to Identify Outliers

### 1️⃣ **Profiling the Data**  
Run `ehdn_new_calls.sh` to generate profiling data using EHDN.

### 2️⃣ **Performing Quality Control**  
Execute `QC.sh` to identify and remove low-quality samples from the **manifest file**.

### 3️⃣ **Merging Data & Running DBSCAN**  
Run `Combine_RUN_DBSCAN.sh` to:
- Aggregate sample data
- Merge and compare regions
- Apply **DBSCAN** for clustering & outlier detection

### 4️⃣ **Annotating Data**  
Run:
- `Annotate_prep_ANNOVAR.r` to prepare the data for annotation
- `annotate_annovar.sh` to annotate using **ANNOVAR**

### 5️⃣ **Cleaning & Preparing Data for Further Analysis**  
Run `prepare_dataframes.r` to refine the dataset for downstream analysis.

### 6️⃣ **Performing Correlation Analysis**  
Run `Correlation_analysis.r` to investigate associations between genomic features and TREs.

### 7️⃣ **Conducting Proximity Analysis**  
Run `proximity_analysis.r` to assess spatial relationships between regulatory elements and TREs.

### 8️⃣ **Executing Statistical Analysis**  
Run `Statistics.r` to perform **Fisher, Kruskal, and Wilcoxon** tests.

### 9️⃣ **Generating Gene Lists**  
Run `gene_lists.r` to extract gene names associated with TREs.

---

## 💻 Command Summary

### 🔹 Run Profiling
```sh
./ehdn_new_calls.sh
```

### 🔹 Perform Quality Control
```sh
./QC.sh -i /home/rachele/EHdn/new_calls/ -m /home/rachele/manifest_noUHRNA.tsv -o /home/rachele/EHDN_DBSCAN_correct/Result/QC
```

### 🔹 Merge Data & Run DBSCAN
```sh
./Combine_RUN_DBSCAN.sh -m /home/rachele/EHDN_DBSCAN_correct/Result/QC/filtered_manifest.tsv
./Combine_RUN_DBSCAN.sh -m /home/rachele/manifest_noUHRNA.tsv
```

### 🔹 Annotate Data
```sh
Rscript Annotate_prep_ANNOVAR.r
./annotate_annovar.sh
```

### 🔹 Clean & Prepare Data
```sh
Rscript prepare_dataframes.r
```

### 🔹 Perform Correlation Analysis
```sh
Rscript Correlation_analysis.r
```

### 🔹 Conduct Proximity Analysis
```sh
Rscript proximity_analysis.r
```

### 🔹 Execute Statistical Analysis
```sh
Rscript Statistics.r
```

### 🔹 Generate Gene Lists
```sh
Rscript gene_lists.r
```

---

## 🎉 Final Thoughts
This pipeline streamlines the process of **outlier detection** in genomic datasets using EHDN and DBSCAN. We hope you find it useful! 🚀

📩 Feel free to contribute or reach out with any questions! 😊









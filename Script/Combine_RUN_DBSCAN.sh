#!/bin/bash

# Usage function
usage() {
    echo "Usage: $0 -m manifest_file"
    exit 1
}

# Parse input arguments
while getopts "i:m:o:" opt; do
    case "$opt" in
        m) manifest_file=$OPTARG ;;
        *) usage ;;
    esac
done

# Check if manifest file is provided
if [ -z "$manifest_file" ]; then
    usage
fi

# Check if the manifest file exists
if [ ! -f "$manifest_file" ]; then
    echo "Error: Manifest file not found: $manifest_file"
    exit 1
fi

# Get the current date and time in the desired format
current_date=$(date +"%Y_%m_%d_%H%M")

# Define the output directory inside the Result folder
output_dir="/home/rachele/EHDN_DBSCAN_correct/Result/results_dbscan_${current_date}"

# Create the output directory
mkdir -p "$output_dir"
if [ $? -ne 0 ]; then
    echo "Error: Failed to create output directory: $output_dir"
    exit 1
fi

# Define the output file names
combined_counts="${output_dir}/EHDn_DBSCAN_${current_date}.combinedCounts.json"
output_regions="${output_dir}/EHDn_DBSCAN_${current_date}.combinedCounts.bed"

# Check if the required Python scripts exist
combine_counts_script="/home/rachele/EHDN_DBSCAN_correct/Script/Fazal2020Scripts/EHDn-v0.6.2_HelperScripts/combine_counts.py"
compare_anchored_irrs_script="/home/rachele/EHDN_DBSCAN_correct/Script/Fazal2020Scripts/EHDn-v0.6.2_HelperScripts/compare_anchored_irrs.py"

if [ ! -f "$combine_counts_script" ]; then
    echo "Error: combine_counts.py script not found: $combine_counts_script"
    exit 1
fi

if [ ! -f "$compare_anchored_irrs_script" ]; then
    echo "Error: compare_anchored_irrs.py script not found: $compare_anchored_irrs_script"
    exit 1
fi

# Run the combine_counts.py script with the current date and time
python3 "$combine_counts_script" --manifest "$manifest_file" --combinedCounts "$combined_counts"
if [ $? -ne 0 ]; then
    echo "Error: combine_counts.py failed"
    exit 1
fi

# Run the compare_anchored_irrs.py script with the current date and time
python3 "$compare_anchored_irrs_script" --manifest "$manifest_file" --inputCounts "$combined_counts" --outputRegions "$output_regions" --minCount 2
if [ $? -ne 0 ]; then
    echo "Error: compare_anchored_irrs.py failed"
    exit 1
fi

# Pass parameters to Rscript via command line arguments
Rscript /home/rachele/EHDN_DBSCAN_correct/Script/dbscan_only.R "$output_dir" "$output_regions" "/home/rachele/EHDN_DBSCAN_correct/wgs.list"
if [ $? -ne 0 ]; then
    echo "Error: dbscan_only.R failed"
    exit 1
fi



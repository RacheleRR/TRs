#!/bin/bash

# Usage function
usage() {
    echo "Usage: $0 -i input_file -m manifest_file -o output_dir"
    exit 1
}

# Parse input arguments
while getopts "i:m:o:" opt; do
    case "$opt" in
        i) input_file=$OPTARG ;;
        m) manifest_file=$OPTARG ;;
        o) output_dir=$OPTARG ;;
        *) usage ;;
    esac
done

# Check if input file, manifest file, and output directory are provided
if [ -z "$input_file" ] || [ -z "$manifest_file" ] || [ -z "$output_dir" ]; then
    usage
fi

# Execute the Python script
python3 /home/rachele/EHDN_DBSCAN_correct/Script/EHDN_raw_call_count_and_tsv.py "$input_file" "$output_dir"

# Execute the R script
Rscript /home/rachele/EHDN_DBSCAN_correct/Script/count_calls.R "$output_dir"

# Filter the manifest file based on the outliers using a Python script
python3 /home/rachele/EHDN_DBSCAN_correct/Script/filter_manifest.py "$output_dir/outliers.tsv" "$manifest_file" "$output_dir/filtered_manifest.tsv" "$output_dir/outlier_df.tsv" "$output_dir/merged_df.tsv"
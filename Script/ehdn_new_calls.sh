#!/bin/bash
# HOW to PROFILE on EHDN

# Define the path to your NAS folder
NAS_FOLDER="/mnt/DRAGEN_pipeline_results"

# Define the reference genome
REF="/home/rachele/Documents/Fasta_files/hg38.fa"

# Define the path to the manifest file
MANIFEST_FILE="/home/rachele/sample.txt"

# Iterate through each sample in the manifest file
while IFS= read -r sample; do
    echo "Processing sample: $sample" 
    
    # Construct the folder path for the current sample
    folder="${NAS_FOLDER}/${sample}/dragen/"

    # Get the BAM file in the sample folder that matches the pattern "*.bam"
    BAM=$(find "$folder" -name "*.bam" -not -name "*.repeats.bam")

    # Check if the BAM file exists
    if [ -f "$BAM" ]; then
        # Define the output directory
        OUTPUT_DIR="/home/rachele/EHdn/new_calls/${sample}"
        mkdir -p "$OUTPUT_DIR"

        # Run ExpansionHunterDenovo on the BAM file
        /home/rachele/Downloads/ExpansionHunterDenovo-v0.9.0-linux_x86_64/bin/ExpansionHunterDenovo profile --read "$BAM" --reference "$REF" --min-anchor-mapq 50 --max-irr-mapq 40 --output-prefix "$OUTPUT_DIR"
                                
        echo "Completed sample: $sample"
    else
        echo "BAM file not found for sample $sample."
    fi

    echo "----------------------------------------" 
done < "$MANIFEST_FILE"

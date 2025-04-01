#!/usr/bin/env python3

import os
import argparse
import BTlib

def process_json_file(json_filename, ignore_alt_contigs, chromosomes):
    try:
        # Open the JSON file
        with open(json_filename) as f:
            sample_count = 0
            counting = False

            # Iterate through each line in the file
            for line in f:
                if "RegionsWithIrrAnchors" in line:
                    # Start counting if the line contains "RegionsWithIrrAnchors"
                    counting = True
                elif "RegionsWithIrrs" in line:
                    # Stop counting if the line contains "RegionsWithIrrs"
                    counting = False
                elif "-" in line and counting:
                    # Extract the chromosome from the line
                    chrom = line.split(":")[0].replace('"', "").replace(" ", "")
                    # Increment the count if conditions are met
                    if not ignore_alt_contigs or chrom.replace("chr", "") in chromosomes:
                        sample_count += 1
        return sample_count
    except IOError as e:
        print(f"Error opening or reading file {json_filename}: {e}")
        return 0

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input_folder", type=str, help="Path to the folder containing JSON files")
    parser.add_argument("--ignore-alt-contigs", default=False, action="store_true")
    args = parser.parse_args()

    # Get the list of chromosomes
    chromosomes = BTlib.get_chromosomes()

    # Initialize a dictionary to store sample counts
    sample_counts = {}

    # Process each JSON file in the specified folder
    for filename in os.listdir(args.input_folder):
        if filename.endswith(".str_profile.json"):
            sample_name = os.path.splitext(filename)[0]
            json_filepath = os.path.join(args.input_folder, filename)
            # Count occurrences in the JSON file
            count = process_json_file(json_filepath, args.ignore_alt_contigs, chromosomes)
            sample_counts[sample_name] = count

    # Write the results to a TSV file
    output_filepath = "sample_counts.tsv"
    try:
        with open(output_filepath, "w") as output_file:
            output_file.write("Sample\tCount\n")
            for sample, count in sample_counts.items():
                output_file.write(f"{sample}\t{count}\n")
        print(f"Results written to {output_filepath}")
    except IOError as e:
        print(f"Error writing to file {output_filepath}: {e}")

if __name__ == "__main__":
    main()
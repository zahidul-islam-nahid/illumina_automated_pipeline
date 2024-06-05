#!/bin/bash

# Set the path to the mother folder containing subfolders A to E
mother_folder="/home/nahid/neisseria/last_two_neisseria"

# Function to filter contig sequences longer than 200 base pairs
filter_contigs() {
    input_fasta="$1"
    output_fasta="$2"
    seqkit seq -m 200 "$input_fasta" > "$output_fasta"
}

#creating file for kraken2
kraken_dir="$mother_folder/kraken"
mkdir -p "$kraken_dir"

# Loop through each subfolder and process the paired-end sequences
for folder in "$mother_folder"/*; do
# Check if the item is a directory
    if [[ -d "$folder" ]]; then
            # Check if the directory is one of the excluded directories
        if [[ "$folder" == */kraken* || "$folder" == */contigs* ]]; then
            echo "Skipping directory $(basename "$folder")."
            continue
        fi
    echo "Processing folder $(basename "$folder")..."

    # Clear the contents of the subfolder, except R1 and R2 files
    find "$folder" -mindepth 1 ! -name "$(basename "$folder")_R1.fastq.gz" ! -name "$(basename "$folder")_R2.fastq.gz" -delete
    
    # Trim the sequences using fastp
    echo "Running fastp for $(basename "$folder")..."
    fastp -i "$folder/$(basename "$folder")_R1.fastq.gz" \
          -I "$folder/$(basename "$folder")_R2.fastq.gz" \
          -f 20 -F 20 -t 10 -T 10 -r -l 40 \
          -o "$folder/$(basename "$folder")_trimmed_R1.fastq" \
          -O "$folder/$(basename "$folder")_trimmed_R2.fastq" | tee "$folder/fastp.log"
    
    # Create a text file containing the trimmed file names for use in fastuniq
    echo "$folder/$(basename "$folder")_trimmed_R1.fastq" >> "$folder/fastuniq.txt"
    echo "$folder/$(basename "$folder")_trimmed_R2.fastq" >> "$folder/fastuniq.txt"
    
    # Use fastuniq to merge the trimmed files
    echo "Running fastuniq for $(basename "$folder")..."
    fastuniq -t q -i "$folder/fastuniq.txt" -o "$folder/$(basename "$folder")_final_R1.fastq" -p "$folder/$(basename "$folder")_final_R2.fastq" | tee "$folder/fastuniq.log"
    
    # Remove the trimmed files
    echo "Deleting trimmed files for $(basename "$folder")..."
    rm "$folder/$(basename "$folder")_trimmed_R1.fastq" "$folder/$(basename "$folder")_trimmed_R2.fastq"
    
    # Compress the final FASTQ files
    echo "Compressing final FASTQ files for $(basename "$folder")..."
    pigz "$folder/$(basename "$folder")_final_R1.fastq" "$folder/$(basename "$folder")_final_R2.fastq"
    
    # Use SPAdes to assemble the merged sequences
    echo "Running SPAdes for $(basename "$folder")..."
    spades -1 "$folder/$(basename "$folder")_final_R1.fastq.gz" -2 "$folder/$(basename "$folder")_final_R2.fastq.gz" -t 26 --cov-cutoff auto --careful -o "$folder/$(basename "$folder")" | tee "$folder/spades.log"

    # Filter contig sequences longer than 200 base pairs
    echo "Filtering contig sequences for $(basename "$folder")..."
    contigs_folder="$mother_folder/contigs"
    mkdir -p "$contigs_folder"
    contigs_fasta="$folder/$(basename "$folder")/contigs.fasta"
    filtered_contigs_fasta="$contigs_folder/$(basename "$folder").fasta"
    filter_contigs "$contigs_fasta" "$filtered_contigs_fasta"
    else
    echo echo "Skipping file $(basename "$folder") as it is not a directory."
    fi
done

# Change the working directory to the contigs directory
contigs_folder="$mother_folder/contigs"
cd "$contigs_folder"

# Run the MLST tool on the filtered contig sequences
echo "Running MLST on filtered contig sequences..."
/home/nahid/tools/mlst/bin/mlst --csv *.fasta  > ../mlst.csv

#running kraken2 for the contig files
for sequence_file in ./*.fasta; do
    filename="${sequence_file##*/}"  # Extract only the file name
    filename="${filename%.fasta}"     # Remove the file extension
    kraken2 --db /home/nahid/tools/kraken2/ --report "${kraken_dir}"/"${filename}".report.txt --report-minimizer-data --output "${kraken_dir}"/"${filename}".output.txt "$sequence_file"
    cut -f1-2,8 "${kraken_dir}"/"${filename}".report.txt > "${kraken_dir}"/"${filename}".std.report.txt
    rm "${kraken_dir}"/"${filename}".report.txt
    rm "${kraken_dir}"/"${filename}".output.txt
done

echo "Script execution completed."

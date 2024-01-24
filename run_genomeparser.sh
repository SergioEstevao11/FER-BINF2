#!/bin/bash

FASTQ_FILE="data/ecoli_PB_small.fastq"
FASTQ_BASENAME=$(basename $FASTQ_FILE)

echo "k-mer_size" > alignment_results.csv

for KMER_SIZE in $(seq 15 2 31); do
    # Run the Python script with the current k-mer size
    python3 genomeparser.py "$FASTQ_FILE" $KMER_SIZE

    # Create a valid filename for the output
    OUTPUT_FILE="minimap2_output_kmer_${KMER_SIZE}_${FASTQ_BASENAME}.txt"

    # Run Minimap2 alignment and store the output in a file
    $(./minimap2/minimap2 -a $FASTQ_FILE output.fa > "$OUTPUT_FILE") > $OUTPUT_FILE

    # Append the k-mer size to the results file
    echo "$KMER_SIZE" >> alignment_results.csv
done

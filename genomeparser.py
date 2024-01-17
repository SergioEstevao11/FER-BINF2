from Bio import SeqIO
import numpy as np
import random
from tqdm import tqdm

num_steps = 100  # You can adjust this value
k = 15
nucleotides = ['A', 'C', 'G', 'T']
states = {}

def parse_genome(filename):
    kmer = ""
    #filename = "data/ecoli_PB_small.fastq"
    with open(filename) as handle:
        for record in tqdm(SeqIO.parse(handle, "fastq"), desc="Building states"):
            seq = record.seq
            for i in tqdm(range(len(seq) - k + 1), desc="Processing kmers", leave=False):
                old_kmer = kmer
                kmer = seq[i:i + k]

                if old_kmer == "":  # first kmer has no previous state, skip
                    continue

                if old_kmer not in states:
                    states[old_kmer] = len(states)

                if kmer not in states:
                    states[kmer] = len(states)

def generate_probabilities(filename):
    num_states = len(states)
    probabilities = np.zeros((num_states, len(nucleotides)))

    # Calculate transition probabilities
    with open(filename) as handle:
        for record in tqdm(SeqIO.parse(handle, "fastq"), desc="Calculating probabilities"):
            seq = record.seq
            for i in tqdm(range(len(seq) - k), desc="Processing kmers", leave=False):
                kmer = seq[i:i + k]
                next_nucleotide = seq[i + k]

                fromId = states[kmer]
                nucleotide_index = nucleotides.index(next_nucleotide)

                probabilities[fromId, nucleotide_index] += 1


    # Normalize the probabilities
    probabilities = probabilities / probabilities.sum(axis=1)[:, np.newaxis]

    return probabilities

def generate_genome(num_steps, probabilities):
    while True:
        initial_state = random.choice(list(states.keys()))
        if any(probabilities[states[initial_state], nucleotide_index] > 0
               for nucleotide_index in range(len(nucleotides))):
            break

    current_state = initial_state
    genome = current_state

    for _ in tqdm(range(num_steps - 1), desc="Generating Genome"):
        nucleotide_index = random.choices(
            range(len(nucleotides)),
            weights=[probabilities[states[current_state], nucleotide_index]
                     for nucleotide_index in range(len(nucleotides))]
        )[0]
        next_nucleotide = nucleotides[nucleotide_index]

        genome += next_nucleotide
        current_state = current_state[1:] + next_nucleotide  # Update the kmer with the new nucleotide

    return genome

def export_genome(genome, output_file="output.fa", header="GEN_genome"):
    with open(output_file, 'w') as file:
        file.write(f">{header}\n")

        genome_str = str(genome) 
        for i in range(0, len(genome_str), 60):
            file.write(genome_str[i:i+60] + "\n")


if __name__ == "__main__":
    gen_filename = "data/ecoli_PB_small.fastq"

    parse_genome(gen_filename)
    probabilities = generate_probabilities(gen_filename)
    generated_genome = generate_genome(num_steps, probabilities)

    export_genome(generated_genome)
    print(generated_genome)

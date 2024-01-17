from Bio import SeqIO
import numpy as np
import random
from tqdm import tqdm

num_steps = 100  # You can adjust this value
k = 15
nucleotides = ['A', 'C', 'G', 'T']

states = {}

# Parse Genome
carry = ""
kmer = ""
filename = "data/ecoli_PB_small.fastq"
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


# Initialize transition probabilities array as a sparse matrix
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

print(probabilities)

def generate_genome(num_steps):
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

# Generate a genome with a given number of steps
generated_genome = generate_genome(num_steps)

# Print or use the generated genome as needed
print(generated_genome)
print(f"genome size: {len(generated_genome)}")
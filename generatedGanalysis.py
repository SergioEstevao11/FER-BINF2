from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from collections import Counter
import tqdm

# File paths
input_genome_file = './data/ecoli_ILL_small.fastq'  # Replace with your file path
output_genome_file = 'output.fa'  # Replace with your file path

# Parameters
k = 15  # k-mer size
nucleotides = ['A', 'C', 'G', 'T']

def parse_genome_for_kmers(filename, k):
    kmer_counts = Counter()
    with open(filename, 'r') as file:
        for record in SeqIO.parse(file, "fastq"):
            seq = record.seq
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i + k]
                kmer_counts[kmer] += 1
    return kmer_counts

def calculate_transition_probabilities(kmer_counts, k):
    states = {}
    state_index = 0
    for kmer in kmer_counts.keys():
        if kmer[:-1] not in states:
            states[kmer[:-1]] = state_index
            state_index += 1

    num_states = len(states)
    transition_matrix = np.zeros((num_states, 4))

    for kmer, count in kmer_counts.items():
        if len(kmer) != k:
            continue
        state = states[kmer[:-1]]
        next_nucleotide = kmer[-1]
        nucleotide_index = nucleotides.index(next_nucleotide)
        transition_matrix[state, nucleotide_index] += count

    # Normalizing the rows to sum to 1
    transition_matrix = np.divide(transition_matrix, transition_matrix.sum(axis=1)[:, np.newaxis], where=transition_matrix.sum(axis=1)[:, np.newaxis] != 0)
    return transition_matrix

# Parsing the input genome file for k-mers
kmer_counts = parse_genome_for_kmers(input_genome_file, k)

# Calculating the transition probabilities
transition_matrix = calculate_transition_probabilities(kmer_counts, k)

# Plotting k-mer frequency distribution
sorted_kmers = sorted(kmer_counts.items(), key=lambda x: x[1], reverse=True)
kmer_labels, kmer_frequencies = zip(*sorted_kmers)

plt.figure(figsize=(10, 6))
plt.bar(kmer_labels[:20], kmer_frequencies[:20])  # Top 20 for visibility
plt.xlabel('K-mer')
plt.ylabel('Frequency')
plt.title('Top 20 K-mer Frequency Distribution')
plt.xticks(rotation=45)
plt.show()
plt.savefig("freq_dist.png")


# Plotting transition probability matrix
plt.figure(figsize=(10, 10))
sns.heatmap(transition_matrix, cmap='viridis')
plt.title('Transition Probability Matrix')
plt.xlabel('To K-mer (Index)')
plt.ylabel('From K-mer (Index)')
plt.show()
plt.savefig("trans_matrix.png")

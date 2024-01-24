import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def plot_kmer_frequency_distribution(kmer_counts):
    sorted_kmers = sorted(kmer_counts.items(), key=lambda x: x[1], reverse=True)
    kmer_labels, kmer_frequencies = zip(*sorted_kmers)

    plt.figure(figsize=(10, 6))
    plt.bar(kmer_labels[:20], kmer_frequencies[:20])  # Top 20 for visibility
    plt.xlabel('K-mer')
    plt.ylabel('Frequency')
    plt.title('Top 20 K-mer Frequency Distribution')
    plt.xticks(rotation=45)
    plt.show()

def plot_transition_probability_matrix(probabilities):
    plt.figure(figsize=(10, 10))
    sns.heatmap(probabilities, cmap='viridis')
    plt.title('Transition Probability Matrix')
    plt.xlabel('To K-mer')
    plt.ylabel('From K-mer')
    plt.show()

def plot_alignment_scores(alignment_scores):
    plt.figure(figsize=(10, 6))
    plt.hist(alignment_scores, bins=30)
    plt.xlabel('Alignment Score')
    plt.ylabel('Frequency')
    plt.title('Distribution of Alignment Scores')
    plt.show()

def plot_genome_coverage(coverage_percentages):
    plt.figure(figsize=(10, 6))
    plt.hist(coverage_percentages, bins=30)
    plt.xlabel('Coverage Percentage')
    plt.ylabel('Frequency')
    plt.title('Genome Coverage Distribution')
    plt.show()

def plot_nucleotide_composition_comparison(original_counts, synthetic_counts):
    labels = list(original_counts.keys())
    original_values = list(original_counts.values())
    synthetic_values = list(synthetic_counts.values())

    x = np.arange(len(labels))  # label locations
    width = 0.35  # width of the bars

    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width/2, original_values, width, label='Original')
    rects2 = ax.bar(x + width/2, synthetic_values, width, label='Synthetic')

    ax.set_xlabel('Nucleotides')
    ax.set_ylabel('Counts')
    ax.set_title('Nucleotide Composition in Original and Synthetic Genomes')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    plt.show()

# Sample data (replace with your actual data)
kmer_counts = {'ACTG': 50, 'CTGA': 45, 'TGAC': 40, 'GACT': 35}  # Example k-mer counts
probabilities = np.random.rand(4, 4)  # Example transition probabilities matrix
alignment_scores = np.random.rand(100) * 100  # Example alignment scores
coverage_percentages = np.random.rand(100) * 100  # Example genome coverage percentages
original_counts = {'A': 300, 'T': 250, 'G': 200, 'C': 150}  # Example original genome nucleotide counts
synthetic_counts = {'A': 290, 'T': 260, 'G': 210, 'C': 140}  # Example synthetic genome nucleotide counts

# Plotting
plot_kmer_frequency_distribution(kmer_counts)
plot_transition_probability_matrix(probabilities)
plot_alignment_scores(alignment_scores)
plot_genome_coverage(coverage_percentages)
plot_nucleotide_composition_comparison(original_counts, synthetic_counts)

from Bio import SeqIO
import itertools
import numpy as np


records = []
states = {}
k = 2
nucleotics = ['A', 'C', 'G', 'T']
probabilities = np.zeros((len(states), len(states)))
filename = "data/ecoli_PB_small.fastq"


# Returns the conditional probability of a transition from a state to another
def getProb(fromId, toId):
    return probabilities[fromId][toId] / sum(probabilities[fromId])


# Initialize possible states
counter = 0
for nuc in itertools.product(nucleotics, repeat=k):
    states[''.join(nuc)] = counter
    counter += 1


# Parse Genome
carry = ""
kmer = ""
with open(filename) as handle:
    for record in SeqIO.parse(handle, "fastq"):
        seq = carry + record.seq # append the end of the last line to the new record
        for i in range(len(seq) - k + 1):
            old_kmer = kmer
            kmer = seq[i:i+k]

            if old_kmer == "": # first kmer has no previous state, skip
                continue         

            fromId = states[old_kmer] 
            toId = states[kmer]

            probabilities[fromId][toId] += 1
            
        carry = seq[(len(seq) - k + 1):len(seq)] # Carries last (k - 1) letters to next line


print(probabilities)

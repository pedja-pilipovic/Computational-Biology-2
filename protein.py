import os
from Bio import SeqIO
import numpy as np
from collections import Counter

class Protein():
    def __init__(self, fasta):
        dirpath = os.getcwd()
        sequences = list(SeqIO.parse(dirpath +'/'+ fasta, "fasta"))
        MSA = []
        for seq in sequences:
            MSA.append(list(str(seq.seq)))
        (N, M) = np.shape(MSA)
        MSA = np.array(MSA)

    def MI(self, sequences, i, j):
        N = np.shape(sequences)[0]
        Pi = Counter(s[i] for s in sequences)
        Pj = Counter(s[j] for s in sequences)
        Pij = Counter((s[i], s[j]) for s in sequences)

        return sum(Pij[(x, y)] / N * np.log(Pij[(x, y)] * N / (Pi[x] * Pj[y])) for x, y in Pij)
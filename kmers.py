import data

class Kmers:
    def __init__(self, reads, k=15):
        self.k = k
        self.reads = data.load_obj(reads)

    def k_mer_counting(self):
        """

      @param reads: List of reads
      @param k: Size of the k-mer by default is 15
      @return: A dictionary which contains the subsequences and the number of times
      that this sequence is repeated k-mer
      """
        kmers = {}
        for N, read in enumerate(self.reads):
            read_str = str(read.seq)
            for start in range(len(read_str) - self.k):
                kmer = read_str[start : start + self.k]
                kmers[kmer] = kmers.get(kmer, 0) + 1
        return kmers


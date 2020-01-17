from Bio import SeqIO
import argparse, os, logging, data

class Kmers:
    def __init__(self, reads, k=15):
        dirpath = os.getcwd()
        self.k = k
        self.reads = list(SeqIO.parse(dirpath + '/' + reads, "fasta"))

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

if __name__ == '__main__':

    #  # define the program description
    text = "This tool takes a simple FASTA files and compute the k-mers algorithm."

    # initiate the parser with a description
    parser = argparse.ArgumentParser(description=text)
    parser.add_argument(
        "-i", "--fasta", help="Name of the FASTA file", required=True
    )
    parser.add_argument(
        "-o", "--output", help="Name of the file to store the k-mers", required=False, default = 'kmers_0'
    )
    parser.add_argument("-k", help="K value", required=False, default=15)
    args = parser.parse_args()

    print('Initializing')
    kmer = Kmers(args.fasta, args.k)
    dictionary_kmers =  kmer.k_mer_counting()
    print('Storing the dictionary of k-mers')
    data.save_obj(dictionary_kmers, args.output + '.pkl')
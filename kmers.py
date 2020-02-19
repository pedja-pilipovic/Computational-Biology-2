"""
Tool that takes a FASTA file and computes the k-mers algorithm.
Export the resulting dictionary into a Pickle file for posterior analysis.
"""

import argparse, os, logging, data
from Bio import SeqIO

__author__ = "BATIER Lucas, GONZALEZ JIMENEZ Alvaro, PILIPVOC Predrag"
__copyright__ = "Copyright 2020, The Salmonela Project"
__license__ = "GPL"
__version__ = "1.0.1"
__email__ = "lucas.batier@hotmail.fr, alvaro.gonzalez-jimenez@grenoble-inp.org, predrag.pilipovic@grenoble-inp.org"
__status__ = "Production"


class Kmers:
    def __init__(self, reads, k=15):
        logging.debug("Initialize the k-mers with k=" + str(k))
        dirpath = os.getcwd()
        self.k = k
        self.reads = list(SeqIO.parse(dirpath + "/" + reads, "fasta"))

    def k_mer_counting(self):
        """

      @param reads: List of reads
      @param k: Size of the k-mer by default is 15
      @return: A dictionary which contains the subsequences and the number of times
      that this sequence is repeated k-mer
      """
        kmers = {}
        logging.debug("K-mer counting")
        for N, read in enumerate(self.reads):
            logging.debug("Subsequences in the read" + str(read.seq))
            read_str = str(read.seq)
            for start in range(len(read_str) - self.k):
                kmer = read_str[start : start + self.k]
                kmers[kmer] = kmers.get(kmer, 0) + 1
        logging.debug("Returning the dictinoary of k-mers")
        return kmers


if __name__ == "__main__":

    #  # define the program description
    text = "This tool takes a simple FASTA files and compute the k-mers algorithm."

    # initiate the parser with a description
    parser = argparse.ArgumentParser(description=text)
    parser.add_argument("-i", "--fasta", help="Name of the FASTA file", required=True)
    parser.add_argument(
        "-o",
        "--output",
        help="Name of the file to store the k-mers",
        required=False,
        default="kmers_0",
    )
    parser.add_argument("-k", help="K value", required=False, default=15)
    args = parser.parse_args()

    print("Initializing")
    kmer = Kmers(args.fasta, int(args.k))
    dictionary_kmers = kmer.k_mer_counting()
    print("Storing the dictionary of k-mers")

    logging.debug("Storing the dictionary of k-mers")
    data.save_obj(dictionary_kmers, args.output + ".pkl")

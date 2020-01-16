# include standard modules
import data, logging, assembler, argparse, math
import numpy as np
from collections import Counter
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

class Snp():
    def __init__(self, kmers_0, kmers_1):
        '''

        :param C: Coverage value
        :param kmers_0: Pickle file of k-mers 0
        :param kmers_1: Pickle file of k-mers 1
        '''
        self.kmers = data.load_obj(kmers_0)
        self.kmers_1 = data.load_obj(kmers_1)
        self.hist_kmers = dict(Counter(kmers_0.values()))

    def mean_coverage(self, start_value=5):
        '''

        :param start_value: Boundary for counting the mean values
        :return: the mean coverage
        '''
        C_mean = 0
        n = 0
        for v in self.hist_kmers:
            if v > start_value:
                C_mean += v
                n += 1
        C_mean /= n

    def pmf_poisson(self,l, k):
        '''

        :param l: lambda, parameters of Poisson law
        :param k: Value
        :return: Poisson probability mass function (PMF)
        '''
        return l ** k * np.exp(-l) / np.math.factorial(k)

    def cdf_poisson(self,l, k):
        '''

        :param l: lambda, parameters of Poisson law
        :param k: Value
        :return: Poisson cumulative density function (CDF)
        '''
        return np.exp(-l) * np.sum([l ** i / np.math.factorial(i) for i in range(0, k + 1)])

    def snp_with_cdf(self, cmean, error_rate=0.01):
        '''

        :param cmean: Mean coverage for the corresponding kmer
        :param error_rate: Minimum error rate to add the SNP
        :return: array of SNPs
        '''
        kmers_0_snp = []
        for kmer in self.kmers.keys():
            if kmer not in self.kmers_1.keys() and self.cdf_poisson(cmean, self.kmers[kmer]) > error_rate:
                kmers_0_snp.append(kmer)
        return kmers_0_snp

if __name__ == '__main__':
    # auxiliar variables to store the path for k-mers files and name for SNPs
    kmers_0 = None
    kmers_1 = None
    snp0_name = None
    snp1_name = None
    k = None
    # define the program description
    text = 'This tool takes two simple k-mers files (from the kmers.py) and outputs a list of SNPs.'

    # initiate the parser with a description
    parser = argparse.ArgumentParser(description=text)
    parser.add_argument("-kmers", "--kmers_file", help="Path of the first k-mers file", required=True)
    parser.add_argument("-kmers2", "--kmers_file2", help="Path of the second k-mers file", required=True)
    parser.add_argument("-snp", "--snp_name", help="Name of the first SNPs file", required=False, default='snp1')
    parser.add_argument("-snp2", "--snp_name2", help="Name of the second SNPs file", required=False, default='snp2')
    parser.add_argument("-k", help="K value", required=False, default=15)

    args = parser.parse_args()

    logging.debug('Getting the path of k-mer file')
    kmers_0 = args.kmers

    logging.debug('Getting the path of the second  k-mer file')
    kmers_1 = args.kmers2

    logging.debug('Getting the name to save the SNPs')
    snp0_name = args.snp

    logging.debug('Getting the name to save the second SNPs')
    snp1_name = args.snp2

    logging.debug('Getting the k value')
    k = args.k

    if (kmers_1 and kmers_0):
        logging.debug('Creating the Snp objects')
        snp_0 = Snp(kmers_0=kmers_0, kmers_1=kmers_1)
        snp_1 = Snp(kmers_0=kmers_1, kmers_1=kmers_0)

        logging.debug('Get the SNPs from the first k-mers file')
        kmers_0_snp = snp_0.snp_with_cdf(cmean=snp_0.mean_coverage())
        logging.info('Save the SNPs from the first k-mers')
        data.save_obj(kmers_0_snp,snp0_name)

        logging.debug('Get the SNPs from the second k-mers file')
        kmers_1_snp = snp_1.snp_with_cdf(cmean=snp_1.mean_coverage())
        logging.info('Save the SNPs from the second k-mers')
        data.save_obj(kmers_1_snp, snp1_name)

        logging.info('Printing the list of SNPs')

        for i in range(math.floor(len(kmers_0_snp) / k)):
            print('SNP ' + str(i + 1))
            # Define two sequences to be aligned
            X = assembler.ah(kmers_0_snp[round(len(kmers_0_snp) / math.floor(len(kmers_0_snp) / k) * i):round(
                len(kmers_0_snp) / math.floor(len(kmers_0_snp) / k) * (i + 1))])
            Y = assembler.ah(kmers_1_snp[round(len(kmers_1_snp) / math.floor(len(kmers_1_snp) / k) * i):round(
                len(kmers_1_snp) / math.floor(len(kmers_1_snp) / k) * (i + 1))])

            # Get a list of the global alignments between the two sequences ACGGGT and ACG
            # No parameters. Identical characters have score of 1, else 0.
            # No gap penalties.
            alignments = pairwise2.align.globalms(X, Y, 2, 0, -1, 0)
            # Use format_alignment method to format the alignments in the list
            print(format_alignment(*alignments[0]))

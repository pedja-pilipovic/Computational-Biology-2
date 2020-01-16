# include standard modules
import getopt, sys, data, logging
import numpy as np
from collections import Counter

# include standard modules
import argparse

class Snp():
    def __init__(self, kmers):
        '''

        :param C: Coverage value
        :param kmers: Pickle file of k-mers
        '''
        logging.debug('Initializing SNP and loading the kmers from: ' + str(kmers))
        self.kmers = data.load_obj(kmers)
        self.hist_kmers = dict(Counter(kmers.values()))

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

    def snp_with_cdf(self, error_rate=0.01):
        '''

        :param error_rate: Minimum error rate to add the SNP
        :return: array of SNPs
        '''
        kmers_snp = []
        for kmer in self.kmers.keys():
            if kmer not in self.kmers.keys() and self.cdf_poisson(self.mean_coverage(), self.kmers[kmer]) > error_rate:
                kmers_snp = np.append(kmers_snp, kmer)
        return kmers_snp

if __name__ == '__main__':
    # auxiliar variables to store the path for k-mers files and name for SNPs
    kmers_1 = None
    snp1_name = None

    # define the program description
    text = 'This script will calculate the SNPs from two differents k-mers files'

    # initiate the parser with a description
    parser = argparse.ArgumentParser(description=text)
    parser.add_argument("-kmers", "--kmers_file", help="Path of the k-mers file")
    parser.add_argument("-snp", "--snp_name", help="Name to store the SNPs")

    args = parser.parse_args()
    try:
        if args.kmers:
            logging.info('Getting the path of k-mer file')
            kmers_1 = args.kmers
        if args.snp:
            logging.info('Getting the name to save the SNPs')
            snp1_name = args.snp

        if (kmers_1 and snp1_name):
            logging.debug('Creating the Snp objects')
            snp_1 = Snp(kmers_1)

            logging.debug('Get the SNPs from the first k-mers file')
            kmers_1_snps = snp_1.snp_with_cdf()
            logging.info('Save the SNPs from the first k-mers')
            data.save_obj(kmers_1_snps, snp1_name)
        else:
            logging.error('Use -h to see how to use the script')
            sys.exit(2)

    except:
        logging.error('Use -h to see how to use the script')

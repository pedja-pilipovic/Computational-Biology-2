# include standard modules
import getopt, sys
import data
import numpy as np
from collections import Counter

class Snp():
    def __init__(self, C, kmers):
        '''

        :param C: Coverage value
        :param kmers: Pickle file of k-mers
        '''
        self.C = C
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


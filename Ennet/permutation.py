#!/usr/bin/env python

'''
Purpose: permutation test

Functions included:
1. count_enh_snps: get snp count of each enhancer


'''

from multiprocessing import cpu_count, Pool
import copy
import networkx as nx
import numpy as np
from scipy.stats import norm

import escore

import logging    # logging

# logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def mean(data):
    """Return the sample arithmetic mean of data.
    http://stackoverflow.com/a/27758326/632242
    """
    n = len(data)
    return sum(data)/float(n)

def _ss(data):
    """Return sum of square deviations of sequence data."""
    c = mean(data)
    ss = sum((x-c)**2 for x in data)
    return ss

def pstdev(data):
    """Calculates the population standard deviation."""
    n = len(data)
    ss = _ss(data)
    pvar = ss/n # the population variance
    return pvar**0.5


def multiple_testing_correction(pvalues, correction_type="Benjamini-Hochberg"):
    """
    Copyright 2017 Francisco Pina Martins <f.pinamartins@gmail.com>

    Taken from https://stackoverflow.com/a/21739593/3091595, remove numpy dependence

    @parameter pvalues - a list of pvalues
    @parameter correction_type - pvalue correction method

    @return qvalues - a list of qvalues
    """
    n = len(pvalues)
    qvalues = [0]*n
    if correction_type == "Bonferroni":
        qvalues = n * pvalues

    elif correction_type == "Bonferroni-Holm":
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        for rank, vals in enumerate(values):
            pvalue, i = vals
            qvalues[i] = (n-rank) * pvalue

    elif correction_type == "Benjamini-Hochberg":
        values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
        values.sort()
        values.reverse()
        new_values = []
        for i, vals in enumerate(values):
            rank = n - i
            pvalue, index = vals
            new_values.append((n/rank) * pvalue)
        for i in xrange(0, int(n)-1):
            if new_values[i] < new_values[i+1]:
                new_values[i+1] = new_values[i]
        for i, vals in enumerate(values):
            pvalue, index = vals
            qvalues[index] = new_values[i]

    return qvalues


def multi_test_pvalues(G):

    '''
    adjust p value

    @parameter G - a graph

    @return G - a graph
    '''

    # get gene pvalues dict
    gene_emp_p = get_value_from_graph(G, 'gene', 'emp_p')

    genes = gene_emp_p.keys()
    emp_p = gene_emp_p.values()

    qvalues = multiple_testing_correction(emp_p)  # multi-test

    gene_qvalue = {genes[i]:qvalues[i] for i in range(len(genes))}

    # put into the graph
    G = put_value_into_graph(gene_qvalue, G, 'gene', 'emp_q')

    return G


def random_net(G):
    '''
    Random the network but keep the original node degree roughly

    @parameter G - graph

    @return GG - random graph that only contains genes
    '''

    H = escore.get_subnetwork(G, 'gene')

    sequence = [d for (_, d) in H.degree]
    GG = nx.configuration_model(sequence)
    GG = nx.Graph(GG)
    GG.remove_edges_from(GG.selfloop_edges())

    # re-label
    mapping = {i:j for i,j in enumerate(H.nodes)}
    GG = nx.relabel_nodes(GG, mapping, copy=False)

    for node in GG.nodes:
        GG.nodes[node]['type'] = 'gene'

    return GG


def permutation_helper(G):
    '''
    helper function for permutation
    '''
    r = G.graph['r']
    p_0 = escore.get_value_from_graph(G, 'gene', 'p_0')
    p_n = escore.get_value_from_graph(G, 'gene', 'p_n')

    GG = random_net(G)
    GG = escore.put_value_into_graph(p_0, GG, 'gene', 'p_0')
    GG.graph['r'] = r
    GG = escore.stationary_p(GG)
    random_p_n = escore.get_value_from_graph(GG, 'gene', 'p_n')

    return random_p_n


def permutation(G, permutation_times, threads):
    '''
    permutation test: random network 500 times to generate emperical distribution of final score

    @parameter G - graph
    @parameter permutation_times - permutation times
    @parameter threads - number of CPUs used to permutate the nework

    @return G - graph
    '''

    if not threads:
        threads = cpu_count()

    logger.info('Start network permutation, using %s threads.' % threads)

    pool = Pool(threads)
    results = []

    # permutation the network to generate emperical p_n for each gene
    for i in range(permutation_times):
        results.append(pool.apply_async(permutation_helper, args=(G,)))

    pool.close()
    pool.join()

    emp_p_n = {}

    for i in results:
        random_p_n = i.get()    # ApplyResult, use get to access the value
        for gene in random_p_n:
            if gene in emp_p_n:
                emp_p_n[gene].append(random_p_n[gene])
            else:
                emp_p_n[gene] = [random_p_n[gene]]

    logger.info('Network permutation done.')

    # test each gene using norm distribution
    for gene in emp_p_n:
        tmp_list = emp_p_n[gene]
        mu, std = mean(tmp_list), pstdev(tmp_list)
        rv = norm(mu, std)
        p = rv.sf(p_n[gene])

        G.nodes[gene]['emp_p'] = p

    logger.info('Computing emperical p-values done.')

    # multi-test adjusting
    G = multi_test_pvalues(G)
    logger.info('Computing adjusted p-values done.')

    return G


def test():
    pass


if __name__ == "__main__":

    test()





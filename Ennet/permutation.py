#!/usr/bin/env python

'''
Purpose: permutation test

Functions included:
1. count_enh_snps: get snp count of each enhancer


'''

import sys
from multiprocessing import cpu_count, Pool
import copy
import networkx as nx
from scipy.stats import norm
import numpy as np
import logging
from numba import jit

import escore
import preprocess

# logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def random_net(G=None):
    '''
    Random the network but keep the original node degree roughly

    @parameter G - graph

    @return GG - random graph that only contains genes
    '''
    sequence = [d for (_, d) in G.degree]
    GG = nx.configuration_model(sequence)
    GG = nx.Graph(GG)
    GG.remove_edges_from(GG.selfloop_edges())

    # re-label
    mapping = {i:j for i,j in enumerate(G.nodes)}
    GG = nx.relabel_nodes(GG, mapping, copy=False)

    for node in GG.nodes:
        GG.nodes[node]['type'] = 'gene'

    return GG


@jit
def get_emp_p(G=None, results=None):
    '''
    get emp_p from permutation test

    @parameter G - graph
    @parameter results - ApplyResult of pool.apply_async

    @return G - graph
    '''
    emp_p_n = dict()
    tmp_list = list()

    for i in results:
        random_p_n = i.get()    # ApplyResult, use get to access the value
        tmp_list.extend(random_p_n)

    tmp_list = np.array(tmp_list)
    total_scores = len(tmp_list)

    # compute p-value of each gene by counting the number of values large than it
    for gene in G:
        G.nodes[gene]['emp_p'] = np.sum(tmp_list>G.nodes[gene]['p_n'])/float(total_scores)

    return G

@jit
def random_stationary_p(A, p_0, r):
    W = escore.normalize_cols(A)
    return list(r*np.dot(np.linalg.inv(np.eye(*np.shape(W))-(1-r)*W), p_0))

def permutation_helper(G=None, p_0=None, r=None):
    '''
    helper function for permutation

    @parameter G - graph
    @parameter p_0 - initial p_0
    @parameter r - restart possibility

    @return random_p_n - p_n on random network
    '''
    GG = random_net(G)
    GG = escore.put_value_into_graph(p_0, GG, 'gene', 'p_0')
    gene_p_0 = escore.get_value_from_graph(GG, 'gene', 'p_0')
    A = nx.to_numpy_array(GG)
    #del GG
    p_0 = np.array(list(gene_p_0.values()))

    p_n = random_stationary_p(A, p_0, r)

    return p_n


def permutation(G=None, permutation_times=None, threads=None):
    '''
    permutation test: random network 500 times to generate emperical distribution of final score

    @parameter G - graph
    @parameter permutation_times - permutation times
    @parameter threads - number of CPUs used to permutate the nework

    @return G - graph
    '''
    if threads is None:
        threads = cpu_count()

    logger.info('Start network permutation, using %s threads.' % threads)

    r = G.graph['r']
    p_0 = escore.get_value_from_graph(G, 'gene', 'p_0')

    pool = Pool(threads)
    results = []
    # permutation the network to generate emperical p_n for each gene
    for i in range(permutation_times):
        results.append(pool.apply_async(permutation_helper, args=(G,p_0,r,)))

    pool.close()
    pool.join()

    logger.info('Network permutation done.')

    G = get_emp_p(G, results)
    logger.info('Computing emperical p-values done.')

    return G


def test():
    network_file = '/home/niuyw/Project/RegulatorySNP_170808/ennet_180821/data/networks/net.txt'
    enhancer_file = '/home/niuyw/Project/RegulatorySNP_170808/ennet_180821/data/interactions/brain.inter.name'
    snp_file = '/home/niuyw/Project/RegulatorySNP_170808/ennet_180821/data/SNPs/brain.rm.hyper'

    G = preprocess.preprocess(network_file, enhancer_file)
    G = escore.escore(snp_file, G, 0.46636)

    G = permutation(G, 10, 8)

    gene_enh_count = escore.get_value_from_graph(G, 'gene', 'enh_num')
    gene_enh_len = escore.get_value_from_graph(G, 'gene', 'enh_len')
    gene_snp_count = escore.get_value_from_graph(G, 'gene', 'snp_count')
    gene_pvalue = escore.get_value_from_graph(G, 'gene', 'pvalue')
    gene_p_0 = escore.get_value_from_graph(G, 'gene', 'p_0')
    gene_p_n = escore.get_value_from_graph(G, 'gene', 'p_n')
    gene_emp_p = escore.get_value_from_graph(G, 'gene', 'emp_p')

    fout = open('permutation_test.txt', 'w')

    for gene in gene_enh_count:
        fout.write('\t'.join([gene, str(gene_enh_count[gene]), str(gene_enh_len[gene]),
                              str(gene_snp_count[gene]), str(gene_pvalue[gene]),
                              str(gene_p_0[gene]), str(gene_p_n[gene]),
                              str(gene_emp_p[gene])]) + '\n')
    fout.close()


if __name__ == "__main__":

    try:
        test()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(0)






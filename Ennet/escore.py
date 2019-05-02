#!/usr/bin/env python

'''
Purpose: calculate escore

Functions included:
1. count_enh_snps: get snp count of each enhancer


'''

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import sys
import networkx as nx
import numpy as np
from scipy.optimize import ridder
from scipy.stats import poisson
from sklearn.preprocessing import normalize
from math import log10

import logging    # logging

# logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def count_enh_snps(G):
    '''
    get snp count of enhancers
    '''

    snp_pos = {contig:{} for contig in G.graph['contigs']}  # a dict of each contig

    with open(G.graph['snp_file'], 'r') as fin:
        for line in fin:
            line = line.strip()
            if line.startswith('#'):
                continue
            line = line.split('\t')
            snp_chr, snp_locus = line[0], line[1]

            if snp_chr not in snp_pos:
                continue

            if snp_locus in snp_pos[snp_chr]:
                snp_pos[snp_chr][snp_locus] += 1
            else:
                snp_pos[snp_chr][snp_locus] = 1

    for node in G:
        if G.nodes[node]['type'] == 'enhancer':
            enh_chr, enh_left, enh_right = node.split('@')
            enh_left, enh_right = int(enh_left), int(enh_right)

            count = 0
            i = enh_left
            while i <= enh_right:
                if not (str(i) in snp_pos[enh_chr]):
                    i += 1
                else:
                    count += snp_pos[enh_chr][str(i)]
                    i += 1
            G.nodes[node]['snp_count'] = count

    return G


def count_gene_snps(G):
    '''
    get snp count of genes
    '''

    # enhancer snp count dict
    enh_snp_count = get_value_from_graph(G, 'enhancer', 'snp_count')

    for x in G:
        if G.nodes[x]['type'] == 'gene':

            if G.nodes[x]['enh_num'] == 0:
                continue

            snp_count = 0
            for enh in G.nodes[x]['enhs']:
                snp_count += enh_snp_count[enh]

            G.nodes[x]['snp_count'] = snp_count

    return G


def get_gene_poisson_p(G):
    '''
    update the gene_pvalue dict

    @parameter gene_snp_count - a gene snp-count dict of the graph
    @parameter G - a graph which contain enhancers and genes

    @return G - a graph
    '''
    # enh snp count dict
    enh_snp_count= get_value_from_graph(G, 'enhancer', 'snp_count')
    G.graph['snp_count'] = sum(enh_snp_count.values())

    # gene snp count dict
    gene_snp_count = get_value_from_graph(G, 'gene', 'snp_count')

    # mutation rate
    mut_rate = G.graph['snp_count']/float(G.graph['enh_len'])

    for x in gene_snp_count:
        enh_len = G.nodes[x]['enh_len']
        snp_count = gene_snp_count[x]

        if snp_count == 0:
            p = 1.0
        else:
            rv = poisson(enh_len*mut_rate) # poisson model pvalue
            p = rv.sf(snp_count-1)   # p-value

        G.nodes[x]['pvalue'] = p

    gene_pvalue = get_value_from_graph(G, 'gene', 'pvalue')
    mini_detected_p = 1.0
    for gene in gene_pvalue:
        if gene_pvalue[gene] == 0:
            continue
        if gene_pvalue[gene] < mini_detected_p:
            mini_detected_p = gene_pvalue[gene]

    mini_detected_p = -log10(mini_detected_p)
    for gene in gene_pvalue:
        if gene_pvalue[gene] == 0:
            G.nodes[gene]['log10_pvalue'] = mini_detected_p
        else:
            G.nodes[gene]['log10_pvalue'] = -log10(gene_pvalue[gene])

    return G


def diffusion_matrix(A, r):
    '''
    Perform the RWR process

    @parameter A - adjacency matrix of a graph
    @parameter r - restart possibility

    @ return diffusion matrix
    '''
    W = normalize(A, norm='l1', axis=0) # column normalized
    return r*np.linalg.inv(np.eye(*np.shape(W))-(1-r)*W)


def difference(A, r):
    '''
    Find difference between fraction of distribution on neighbors and non-neighbors

    @parameter A - adjacency matrix of a graph
    @parameter r - restart possibility

    @return - difference between fraction of distribution on neighbors and non-neighbors
    '''
    F = diffusion_matrix(A, r)
    np.fill_diagonal(F, 0)

    n = np.sum(F[np.where(A>=1)])    # non-neighbors
    s = np.sum(F[np.where(A<1)])    # neighbors
    return n-s


def choose_r(G):
    '''
    Find value of r that sets difference to zero between fraction of distribution on neighbors and non-neighbors to zero.
    Derived from Hotnet2.

    @parameter G - a graph

    @return G - a graph
    '''
    A = nx.to_numpy_array(G)
    r = ridder(lambda r: difference(A, r), a=0.01, b=0.99, xtol=0.001)

    G.graph['r'] = r
    return G


def plot_root_finding(G):
    '''
    plot the process of root finding

    @parameter G - a graph
    '''

    A = nx.to_numpy_array(G)

    f = lambda x: difference(A, x)
    x = np.linspace(0, 1, 100)
    y = f(x)
    print(x, y)
    plt.plot(x, y)
    plt.axhline(0, color='k')
    plt.xlim(0, 1)

    plt.axis('off')
    plt.savefig("roo_finding.png", dpi=300)
    plt.show()


def get_p_0(G, score_method='log10_pvalue'):
    '''
    get init p0 of RWR

    @parameter G - graph
    @parameter score_method - scoring method of p0: snp_count, pvalue, log10_pvalue

    @return G - graph
    '''

    if score_method == 'snp_count':
        gene_snp_count = get_value_from_graph(G, 'gene', 'snp_count')
        total_value = float(sum(gene_snp_count.values()))

        for gene in gene_snp_count:
            G.nodes[gene]['p_0'] = gene_snp_count[gene]/total_value

    elif score_method == 'pvalue':
        gene_pvalue = get_value_from_graph(G, 'gene', 'pvalue')
        total_value = sum([1-i for i in gene_pvalue.values()])

        for gene in gene_pvalue:
            G.nodes[gene]['p_0'] = (1-gene_pvalue[gene])/total_value

    elif score_method == 'log10_pvalue':
        gene_log10_pvalue = get_value_from_graph(G, 'gene', 'log10_pvalue')
        total_value = float(sum(gene_log10_pvalue.values()))

        for gene in gene_log10_pvalue:
            G.nodes[gene]['p_0'] = gene_log10_pvalue[gene]/total_value

    return G


def stationary_p(G):
    '''
    Calculate p when RWR reaches a stationary distribution

    @parameter G - graph

    @return G - graph
    '''
    A = nx.to_numpy_array(G)
    W = normalize(A, norm='l1', axis=0)
    r = G.graph['r']

    gene_p_0 = get_value_from_graph(G, 'gene', 'p_0')
    p_0 = np.array(list(gene_p_0.values()))

    p_n = list(r*np.dot(np.linalg.inv(np.eye(*np.shape(W))-(1-r)*W), p_0))

    for i, gene in enumerate(list(gene_p_0.keys())):
        G.nodes[gene]['p_n'] = p_n[i]

    return G


def put_value_into_graph(value_dict, G, node_type, value_type):
    '''
    put values into a graph

    @parameter value_dict - a key-value dict
    @parameter G - a graph
    @parameter node_type - the node type
    @parameter value_type - the value type

    @return - a a graph
    '''

    H = G.copy()
    for node in value_dict:
        if H.nodes[node]['type'] == node_type:
            H.nodes[node][value_type] = value_dict[node]

    return H


def get_value_from_graph(G, node_type, value_type):
    '''
    get values from a graph

    @parameter G - a graph
    @parameter node_type - the node type
    @parameter value_type - the value type

    @return - a dict of key-value
    '''
    value_dict = {}

    for node in G:
        if G.nodes[node]['type'] == node_type:
            value_dict[node] = G.nodes[node][value_type]

    return value_dict


def escore(snp_file, G, r=None):
    '''
    escore

    @snp_file - a bed file for SNP locations
    @G - a graph

    @return - a graph
    '''

    # count the SNPs of enhancers
    G.graph['snp_file'] = snp_file
    G = count_enh_snps(G)
    logger.info('Counting SNPs in enhancers done.')

    # count the SNPs of genes
    G = count_gene_snps(G)
    logger.info('Counting SNPs of genes done.')

    # Poisson test to get p value
    G = get_gene_poisson_p(G)
    logger.info('Computing Poisson p-values done.')

    # only keep genes
    gene_snp_count = get_value_from_graph(G, 'gene', 'snp_count')
    G = G.subgraph(list(gene_snp_count.keys()))

    # choose r
    if r:
        G.graph['r'] = r
        logger.info('Using restart possibility provided: %s' % r)
    else:
        logger.info('Restart possibility r not provided. Will be computed from the network.')
        G = choose_r(G)
        logger.info('Choosing restart possibility done: r is: %s' % G.graph['r'])

    #plot_root_finding(G)

    # get p_0
    G = get_p_0(G, 'log10_pvalue')

    # get stationary p
    G = stationary_p(G)
    logger.info('Computing stationary possibilities of nodes done.')

    return G


def test():
    import preprocess

    network_file = '/home/niuyw/Project/RegulatorySNP_170808/ennet_180821/data/networks/net.txt'
    enhancer_file = '/home/niuyw/Project/RegulatorySNP_170808/ennet_180821/data/interactions/brain.inter.name'
    snp_file = '/home/niuyw/Project/RegulatorySNP_170808/ennet_180821/data/SNPs/brain.rm.hyper'

    G = preprocess.preprocess(network_file, enhancer_file)

    G = escore(snp_file, G)

    gene_enh_count = get_value_from_graph(G, 'gene', 'enh_num')
    gene_enh_len = get_value_from_graph(G, 'gene', 'enh_len')
    gene_snp_count = get_value_from_graph(G, 'gene', 'snp_count')
    gene_pvalue = get_value_from_graph(G, 'gene', 'pvalue')
    gene_p_0 = get_value_from_graph(G, 'gene', 'p_0')
    gene_p_n = get_value_from_graph(G, 'gene', 'p_n')

    fout = open('escore_test.txt', 'w')

    for gene in gene_enh_count:
        fout.write('\t'.join([gene, str(gene_enh_count[gene]), str(gene_enh_len[gene]),
                              str(gene_snp_count[gene]), str(gene_pvalue[gene]),
                              str(gene_p_0[gene]), str(gene_p_n[gene])]) + '\n')
    fout.close()


if __name__ == "__main__":

    try:
        test()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(0)




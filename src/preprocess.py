#!/usr/bin/env python

'''
Purpose: preprocess the inputs

Functions included:
1. init_graph: read the network file, init the graph.
2. get_gene_enh_pair: read the interaction file, put the enhancer and gene into the graph.
3. preprocess: do the jobs above.

'''


import os
import platform
import networkx as nx

import logging    # logging

# logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def init_graph(network_file):

    G = nx.Graph(network_file = os.path.abspath(network_file))
    G.graph['python_version'] = platform.python_version()

    with open(network_file, 'r') as fin:
        for line in fin:
            line = line.strip()

            if line.startswith('#'):    # skip lines starting with #
                continue

            gene1, gene2 = line.split()[:2]
            G.add_edge(gene1, gene2, type='gg') # gg for gene-gene interaction

            # add attributes
            G.nodes[gene1]['type'] = 'gene'; G.nodes[gene2]['type'] = 'gene'  # node type as gene
            G.nodes[gene1]['snp_count'] = 0; G.nodes[gene2]['snp_count'] = 0  # raw snp count
            G.nodes[gene1]['enhs'] = set(); G.nodes[gene2]['enhs'] = set()  # store enhancers
            G.nodes[gene1]['enh_len'] = 0; G.nodes[gene2]['enh_len'] = 0  # total enhancer length

    # remove self-conneted nodes
    G.remove_edges_from(G.selfloop_edges())

    return G


def get_gene_enh_pair(enhancer_file, G):
    '''
    get enhancer gene pairs
    '''

    G.graph['enhancer_file'] = os.path.abspath(enhancer_file)
    G.graph['enh_len'] = 0  # total enhancer length
    G.graph['enh_num'] = 0  # total enhancer number
    G.graph['snp_count'] = 0  # total SNP count

    with open(enhancer_file, 'r') as fin:
        for line in fin:
            line = line.strip()
            if line.startswith('#'):
                continue

            line = line.split('\t')

            enhancer, enh_target = '@'.join(line[:3]), line[3]
            enh_len = int(line[2]) - int(line[1])

            G.graph['enh_len'] += enh_len   # update the total enhancer length
            G.graph['enh_num'] += 1 # update the total enhancer number

            G.add_edge(enhancer, enh_target, type='eg') # eg for enhancer-gene

            G.nodes[enhancer]['type'] = 'enhancer'; G.nodes[enh_target]['type'] = 'gene'   # node type as gene
            G.nodes[enhancer]['snp_count'] = 0; G.nodes[enh_target]['snp_count'] = 0   # raw snp count
            G.nodes[enhancer]['enh_len'] = enh_len  # enhancer length

            if 'enhs' not in G.nodes[enh_target]:
                G.nodes[enh_target]['enhs'] = set([enhancer])
            else:
                G.nodes[enh_target]['enhs'].add(enhancer)

            if 'enh_len' not in G.nodes[enh_target]:
                G.nodes[enh_target]['enh_len'] = enh_len
            else:
                G.nodes[enh_target]['enh_len'] += enh_len

    return G


def preprocess(network_file, enhancer_file):

    G = init_graph(network_file)
    x1, y1 = G.number_of_nodes(), G.number_of_edges()

    logger.info('Network initiation done: %s genes and %s links included.' %(x1, y1))

    G = get_gene_enh_pair(enhancer_file, G)
    x2, y2 = G.number_of_nodes(), G.number_of_edges()

    enh_num = 0
    for i in G:
        if G.nodes[i]['type']=='enhancer':
            enh_num += 1

    new_gene_num = 0
    for i in G:
        if G.nodes[i]['type']=='gene':
            if G.nodes[i]['enh_len']:
                new_gene_num += 1

    logger.info('Reading enhancer-gene interactions done: %s enhancers, %s genes and %s interactions included. And %s genes are not in the network.' %(enh_num, new_gene_num, y2-y1, x2-x1-enh_num))

    return G


def test():
    network_file = '/home/niuyw/Project/RegulatorySNP_170808/ennet_180821/data/networks/net.txt'
    enhancer_file = '/home/niuyw/Project/RegulatorySNP_170808/ennet_180821/data/interactions/breast.ep'

    G = preprocess(network_file, enhancer_file)


if __name__ == "__main__":

    test()




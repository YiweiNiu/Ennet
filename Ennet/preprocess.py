#!/usr/bin/env python

'''
Purpose: preprocess the inputs

Functions included:
1. init_graph: read the network file, init the graph.
2. get_gene_enh_pair: read the interaction file, put the enhancer and gene into the graph.
3. preprocess: do the jobs above.

'''

import os
import sys
import platform
import networkx as nx
import logging    # logging

import escore

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
            G.nodes[gene1]['enh_num'] = 0; G.nodes[gene2]['enh_num'] = 0  # total enhancer length

    # remove self-conneted nodes
    G.remove_edges_from(G.selfloop_edges())

    return G


def get_gene_enh_pair(enhancer_file, G):
    '''
    get enhancer gene pairs
    '''

    G.graph['enhancer_file'] = os.path.abspath(enhancer_file)
    G.graph['contigs'] = set()

    with open(enhancer_file, 'r') as fin:
        for line in fin:
            line = line.strip()
            if line.startswith('#'):
                continue

            line = line.split('\t')

            contig, enh_target = line[0], line[3]

            # exclude genes not in the gene-gene network
            if enh_target not in G:
                continue

            G.graph['contigs'].add(contig)  # store all the contigs, used in counting snps

            enhancer = '@'.join(line[:3])
            enh_len = int(line[2]) - int(line[1])

            G.add_edge(enhancer, enh_target, type='eg') # eg for enhancer-gene

            G.nodes[enhancer]['type'] = 'enhancer'
            G.nodes[enhancer]['snp_count'] = 0
            G.nodes[enhancer]['enh_len'] = enh_len

            G.nodes[enh_target]['enhs'].add(enhancer)  # store enhs

    return G


def preprocess(network_file, enhancer_file):

    G = init_graph(network_file)
    x1, y1 = G.number_of_nodes(), G.number_of_edges()

    logger.info('Network initiation done: %s genes and %s links included.' %(x1, y1))

    G = get_gene_enh_pair(enhancer_file, G)
    y2 = G.number_of_edges()

    # get enh lens
    enh_lens = escore.get_value_from_graph(G, 'enhancer', 'enh_len')
    G.graph['enh_len'] = sum(enh_lens.values())  # total enhancer length
    G.graph['enh_num'] = len(enh_lens)  # total enhancer number

    # update genes' enh_len and enh_num
    gene_with_enhs = 0
    for x in G:
        if G.nodes[x]['type'] == 'gene':
            tmp = len(G.nodes[x]['enhs'])
            G.nodes[x]['enh_num'] = tmp

            if tmp != 0:
                gene_with_enhs += 1

                for enh in G.nodes[x]['enhs']:
                    G.nodes[x]['enh_len'] += G.nodes[enh]['enh_len']

    logger.info('Reading enhancer-gene interactions done: %s enhancers and %s interactions included. %s genes have at least one enhancer.' %(G.graph['enh_num'], y2-y1, gene_with_enhs))

    return G


def test():
    network_file = '/home/niuyw/Project/RegulatorySNP_170808/ennet_180821/data/networks/net.txt'
    enhancer_file = '/home/niuyw/Project/RegulatorySNP_170808/ennet_180821/data/interactions/brain.inter.name'

    G = preprocess(network_file, enhancer_file)


if __name__ == "__main__":

    try:
        test()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(0)





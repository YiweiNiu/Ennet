#!/usr/bin/env python
import numpy as np
import networkx as nx
import os
import platform
from sklearn.preprocessing import normalize
from numba import jit
from scipy.optimize import ridder


def init_graph(network_file=None):

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


def normalize_cols(A=None):
    '''
    normalize matrix by cols
    
    @parameter A - adjacency matrix of a graph

    @return - column-normalized matrix
    '''
    return normalize(A, norm='l1', axis=0)


@jit
def diffusion_matrix(A=None, r=None):
    '''
    Perform the RWR process

    @parameter A - adjacency matrix of a graph
    @parameter r - restart possibility

    @ return diffusion matrix
    '''
    W = normalize_cols(A) # column normalized
    return r*np.linalg.inv(np.eye(*np.shape(W))-(1-r)*W)


@jit
def difference(A=None, r=None):
    '''
    Find difference between fraction of distribution on neighbors and non-neighbors

    @parameter A - adjacency matrix of a graph
    @parameter r - restart possibility

    @return - difference between fraction of distribution on neighbors and non-neighbors
    '''
    F = diffusion_matrix(A, r)
    t = np.trace(F)

    np.fill_diagonal(F, 0)

    n = np.sum(F[np.where(A>=1)])    # non-neighbors
    s = np.sum(F[np.where(A<1)])    # neighbors
    #return n - s
    return 0.5*t - n - s


# numba does not support ridder function
def choose_r(G=None):
    '''
    Find value of r that sets difference to zero between fraction of distribution on neighbors and non-neighbors to zero.
    Derived from Hotnet2.

    @parameter G - a graph

    @return G - a graph
    '''
    A = nx.to_numpy_array(G)
    r = ridder(lambda r: difference(A, r), a=0.01, b=0.99, xtol=0.001)

    return r

if __name__ == '__main__':
    G = init_graph("/home/niuyw/Project/RegulatorySNP_170808/ennet_190429/data/networks/net.txt")
    r = choose_r(G)
    print("r="+str(r))

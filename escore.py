#!/usr/bin/env python

'''
Purpose: given an enhancer-gene interactions and a somatic mutation file, return enhancer_snp_count and gene_snp_count

Functions included:
1. count_enh_snps: get snp count of each enhancer


'''

import networkx as nx
import scipy.stats as st
import math
from sklearn.preprocessing import normalize
from scipy.linalg import inv
import scipy as sp
import numpy as np
from math import log10


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


def count_enh_snps(snp_file, G):
    '''
    get snp count of enhancers
    '''

    snp_pos = {"chr1": {}, "chr2": {}, "chr3": {}, "chr4": {}, "chr5": {},
               "chr6": {}, "chr7": {}, "chr8": {}, "chr9": {}, "chr10": {},
               "chr11": {}, "chr12": {}, "chr13": {}, "chr14": {}, "chr15": {},
               "chr16": {}, "chr17": {}, "chr18": {}, "chr19": {}, "chr20": {},
               "chr21": {}, "chr22": {}, "chrX": {}, "chrY": {}}

    with open(snp_file, 'r') as fin:
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

            snp_count = 0
            x_neighbors = G.neighbors(x)
            G.nodes[x]['enhcount'] = 0

            # count total snps of a gene
            for y in x_neighbors:
                if G.nodes[y]['type'] == 'enhancer':
                    G.nodes[x]['enhcount'] += 1
                    snp_count += enh_snp_count[y]

            G.nodes[x]['snp_count'] = snp_count

    return G

def get_weight_and_firstp(G):
    H=G.copy()
    #seed_nodes=[]
    whole_weight=0
    node_need_remove=[]
    for node in H:
        if H.nodes[node]['type'] == 'enhancer':
            node_need_remove.append(node)
        #elif 'snp_count' in H.nodes[node]:
            #if H.nodes[node]['snp_count']>0:
                #seed_nodes.append(node)
                #whole_weight += H.nodes[node]['snp_count']
        elif H.degree()[node] < 1:
            node_need_remove.append(node)
        else:
            whole_weight += (1 - H.nodes[node]['pvalue'])
    H.remove_nodes_from(node_need_remove)

    init_w=nx.to_numpy_matrix(H)
    norm_w=normalize(init_w,norm='l1',axis=0) #normalize with each column(sum=1)
    
    p_0=[0] * H.number_of_nodes()
    for node in H:
        node_index=list(H.nodes).index(node)
        p_0[node_index]=(1 - H.nodes[node]['pvalue'])/whole_weight
    p_0=np.array(p_0)
    return H, init_w, norm_w, p_0

def similarity_matrix(A, r):
    '''
    Perform the random walk with restart process in Ennet
    '''
    W=normalize(A,norm='l1',axis=0) #normalize with each column(sum=1)
    return r*inv(np.eye(*np.shape(W))-(1-r)*W)

def difference(A, r):
    '''
    Find difference between fraction of distribution on neighbors and non-neighbors.
    '''
    P = similarity_matrix(A, r)
    np.fill_diagonal(P, 0)

    n = np.sum(P[np.where(A>=1)])
    s = np.sum(P[np.where(A<1)])
    return n-s
    

def choose_r(A):
    '''
    Find value of r that sets difference to zero between fraction of distribution on neighbors and non-neighbors to zero. Derived from Hotnet2.
    '''
    return sp.optimize.ridder(lambda r: difference(A, r), a=0.01, b=0.99, xtol=0.001)

def iterate_p(p,weight_matrix,r,p_0):
    p_new=np.dot(weight_matrix,p)*(1-r)+p_0*r
    diff_norm=np.linalg.norm(np.subtract(p_new, p), 1)
    if diff_norm <= 0.000001:
        return p
    else:
        #print(p)
        return iterate_p(p_new,weight_matrix,r,p_0)

def stationary_p(weight_matrix,r,p_0):
    '''
    Calculate p when it reaches a stationary distribution
    '''
    return r*np.dot(inv(np.eye(*np.shape(weight_matrix))-(1-r)*weight_matrix),p_0)


def get_gene_poisson_p(G):
    '''
    update the gene_pvalue dict's list

    @parameter gene_snp_count - a gene snp-count dict of the graph
    @parameter G - a graph which contain enhancers and genes

    @return G - a graph
    '''

    # gene snp count dict
    gene_snp_count = get_value_from_graph(G, 'gene', 'snp_count')

    # mutation rate
    mut_rate = sum(gene_snp_count.values())/float(G.graph['enh_len'])

    for x in gene_snp_count:
        enh_len = G.nodes[x]['enh_len']
        snp_count = gene_snp_count[x]

        if snp_count == 0:
            p = 1.0
        else:
            rv = st.poisson(enh_len*mut_rate) # poisson model pvalue
            p = 1 - rv.cdf(snp_count - 1)   # p-value

        G.nodes[x]['pvalue'] = p

    return G


def combine_neighbor_pvalues(G):
    '''
    combine pvalues of neighbor nodes as the final p-value of a node using Fisher's method

    @parameter G - a graph

    @return G - a graph
    '''

    H = G.copy()
    for x in H:

        if G.nodes[x]['type'] == 'gene':

            x_neighbors = [i for i in G.neighbors(x) if G.nodes[i]['type']=='gene']

            if not x_neighbors:     # if one node doesn't have a neighbor gene, it will not be outputed
                H.nodes[x]['pvalue'] = 1.0
                continue

            tmp_pvalues = [G.nodes[x]['pvalue']] # inlcude the node itself

            for y in x_neighbors:
                tmp_pvalues.append(G.nodes[y]['pvalue']) # all neighbors pvalues

            H.nodes[x]['pvalue'] = st.combine_pvalues(tmp_pvalues)[1]   # combine pvalues use Fisher's method

    return H


def multi_test_pvalues(G):

    '''
    combine pvalues of neighbor nodes as the final p-value of a node using Fisher's method

    @parameter G - a graph

    @return G - a graph
    '''

    # get gene pvalues dict
    gene_pvalue = get_value_from_graph(G, 'gene', 'pvalue')

    genes = gene_pvalue.keys()
    pvalues = gene_pvalue.values()

    qvalues = multiple_testing_correction(pvalues)  # multi-test

    x = xrange(len(genes))

    gene_qvalue = {genes[i]:qvalues[i] for i in x}

    # put into the graph
    G = put_value_into_graph(gene_qvalue, G, 'gene', 'qvalue')

    return G


def put_value_into_graph(tmp_dict, G, node_type, value_type):
    '''
    put values into a graph

    @parameter tmp_dict - a value dict
    @parameter G - a graph
    @parameter node_type - the node type
    @parameter value_type - the value type

    @return - a a graph
    '''

    H = G.copy()
    for node in tmp_dict:
        if H.nodes[node]['type'] == node_type:
            H.nodes[node][value_type] = tmp_dict[node]

    return H


def get_value_from_graph(G, node_type, value_type):
    '''
    get values from a graph

    @parameter G - a graph
    @parameter node_type - the node type
    @parameter value_type - the value type

    @return - a dict of gene-pvalues
    '''
    tmp_dict = {}

    for node in G:
        if G.nodes[node]['type'] == node_type:
            tmp_dict[node] = G.nodes[node][value_type]

    return tmp_dict


def escore(snp_file, G):
    '''
    escore
    '''

    G = count_enh_snps(snp_file, G)
    G = count_gene_snps(G)

    G = get_gene_poisson_p(G)

    #G = combine_neighbor_pvalues(G)

    #G = multi_test_pvalues(G)

    G, aj_matrix, weight_matrix, p_0 = get_weight_and_firstp(G)
    r=choose_r(aj_matrix)
    print("Choose r=%s" % (r))

    p_n = stationary_p(weight_matrix,r,p_0)
    #p_n = iterate_p(p_0,weight_matrix,r,p_0)
    
    i = 0
    for node in G:
        G.nodes[node]['staypos'] = p_n[i]
        i = i + 1


    return G


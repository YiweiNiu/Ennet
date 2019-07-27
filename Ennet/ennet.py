#!/usr/bin/env python

'''
Purpose: main function of ennet.

'''

import os
os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=1
os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=1
os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=1
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # export VECLIB_MAXIMUM_THREADS=1
os.environ["NUMEXPR_NUM_THREADS"] = "1" # export NUMEXPR_NUM_THREADS=1

import sys
import argparse
from datetime import datetime
import networkx as nx
from collections import OrderedDict
import logging
import pickle

import preprocess
import escore
import permutation

# logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

positive_genes=["ACTRT1","ASCL1","ATF7IP","BHLHE41","BMF","BMP4","BTN3A2","CCDC171",
                "CCND1","CCNE1","CD274","CD47","CEACAM21","CXCR4","ELL2","EPHB3",
                "ESR1","ETV1","FGF10","FOXE1","GFI1B","GREM1","HOTAIR","IGFBP5",
                "IL1B","LMO1","LUNAR1","MECOM","MGMT","MIR1262","MLH1","MRPS30",
                "MUC5AC","MUC5B","MYB","MYC","PAX5","PCAT1","PCAT19","PDE4B","PSIP1",
                "PTCSC2","PTX3","PVT1","RARA","RASL11A","RFX6","RMND1","SMAD7","SOD2",
                "SOX9","ST18","STC1","TAL1","TERT","TP63","UPK3A"]


def report(G=None, output_prefix=None):
    '''
    '''
    gene_enh_count = escore.get_value_from_graph(G, 'gene', 'enh_num')
    gene_enh_len = escore.get_value_from_graph(G, 'gene', 'enh_len')
    gene_snp_count = escore.get_value_from_graph(G, 'gene', 'snp_count')
    gene_pvalue = escore.get_value_from_graph(G, 'gene', 'pvalue')
    gene_p_0 = escore.get_value_from_graph(G, 'gene', 'p_0')
    gene_p_n = escore.get_value_from_graph(G, 'gene', 'p_n')
    gene_emp_p = escore.get_value_from_graph(G, 'gene', 'emp_p')

    total_enh_len = G.graph['enh_len']
    total_enh_snp_count = G.graph['snp_count']

    gene_emp_p_orig = sorted(gene_emp_p.items(), key=lambda x:x[1])
    gene_emp_p_sorted_dict = OrderedDict()
    for pair in gene_emp_p_orig:
        gene_emp_p_sorted_dict[pair[0]]=pair[1]

    nodeList = open('%s_nodes.txt' %(output_prefix), 'w')

    i = 1
    positive_gene_rank={}
    for node in gene_emp_p_sorted_dict:
        if node in positive_genes:
            positive_gene_rank[node]=i
        # gene name, total enhancer length, total snp count, gene enhancer length
        # gene snp count, gene pvalue, gene_p_0, gene_p_n, gene_emp_p, gene_emp_q
        nodeList.write('\t'.join([node, str(total_enh_len), str(total_enh_snp_count),
                                  str(gene_enh_len[node]), str(gene_snp_count[node]),
                                  str(gene_pvalue[node]), str(gene_p_0[node]),
                                  str(gene_p_n[node]), str(gene_emp_p_sorted_dict[node])]) + '\n')
        i += 1

    nodeList.close()

    positive_gene_detail = open('%s_positive_detail.txt' %(output_prefix), 'w') #positive gene enhancer number and snp count
    for gene in positive_genes:
        if gene not in positive_gene_rank:
            continue
        if gene not in gene_enh_len:
            gene_enh_len[gene]=0
        if gene not in gene_snp_count:
            gene_snp_count[gene]=0
        positive_gene_detail.write('\t'.join([gene, str(gene_enh_count[gene]), str(gene_enh_len[gene]),
                                              str(gene_snp_count[gene]), str(positive_gene_rank[gene])])+'\n')
    positive_gene_detail.close()

    network_object_file = open('%s_network.pickle' %(output_prefix), 'wb')
    pickle.dump(G, network_object_file, 2)
    network_object_file.close()
    logger.info('Writing results done.')


def ennet(args):
    '''
    '''
    network_file = args.network
    enhancer_file = args.enhancer
    snp_file = args.mutation
    # do some args checking here

    # preprocess
    G = preprocess.preprocess(network_file, enhancer_file)

    # escore
    G = escore.escore(snp_file, G, args.r, args.p0_meth)

    # permutation test
    G = permutation.permutation(G, args.permutation, args.threads)

    # report results
    report(G, args.output)

    logger.info('All done. Cheers')


def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, \
                                    usage='\n\npython runEnnet.py [-h] <-n network file> <-e enhancer-gene pairs> <-m snp file> [-r restart possibility] [-o output prefix]\n', \
                                    description='', epilog='Ennet: enhancer-mutated cancer gene prioritizaion with random walk with restart.')

    parser.add_argument('-t', '--threads', type=int, help='Number of CPUs used when permutation. Using all CPUs by default.', metavar='', dest="threads")
    parser.add_argument('-n', '--network', required=True, type=str, help='Path to tab-separated network.', metavar='', dest="network")
    parser.add_argument('-e', '--enhancer', required=True, type=str, help='Path to enhancer-gene pairs.', metavar='', dest="enhancer")
    parser.add_argument('-m', '--mutation', required=True, type=str, help='Path to snp file.', metavar='', dest="mutation")
    parser.add_argument('-r', type=float, help='Restart possibility. Choose between 0 and 1.', metavar='', dest="r")
    parser.add_argument('-p', '--permutation', default=500, type=int, help='Permutation times.', metavar='', dest="permutation")
    parser.add_argument('-o', '--output', default='ennet_res', help='Output prefix.', metavar='', dest="output")
    parser.add_argument('--p0-method', default='log10_pvalue', choices=['log10_pvalue', 'snp_count', 'pvalue'], type=str, help='Method to compute p0.', metavar='', dest='p0_meth')

    if len(sys.argv) < 2:
        argparser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    ennet(args)


if __name__ == '__main__':

    try:
        main()
    except KeyboardInterrupt:
        logger.error("User interrupted me! ;-) Bye!")
        sys.exit(0)


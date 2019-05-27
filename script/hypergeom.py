#!/usr/bin/env python
# ./hypergeom.py result_file true_file

import sys
from scipy import stats
import pandas as pd
import logging


if __name__ == '__main__':

    logging.basicConfig(level = logging.INFO, format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    logging.info('Loading input')
    true_file = open(sys.argv[2], 'r')
    true_gene = []
    for line in true_file:
        line = line.strip()
        true_gene.append(line)
    true_number = len(true_gene)
    result_file = pd.read_csv(sys.argv[1], sep = '\t', header = None)
    find_gene = list(result_file.loc[result_file[8] <= 0.01,].loc[:, 0])
    whole_number = result_file.shape[0]
    find_number = len(find_gene)
    find_true_gene = []
    for gene in find_gene:
        if gene in true_gene:
            find_true_gene.append(gene)
    find_true_number = len(find_true_gene)
    logging.info('Total %d genes, %d true genes. Find %d genes, %d true genes.' % (whole_number, true_number, find_number, find_true_number))
    pvalue = stats.hypergeom.sf(find_true_number, whole_number, find_number, true_number)
    logging.info('The Hypergeometric p value is %f' % (pvalue))
    true_file.close()
    


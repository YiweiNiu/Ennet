#!/usr/bin/env python

import sys
import pickle
import networkx as nx
from Ennet import permutation
from Ennet import ennet

G = pickle.load(open(sys.argv[1], 'rb'))

H = G.copy()

tmp_list = []
for gene in H:
    tmp_list.extend(H.nodes[gene]['emp_p_n'])

total_scores = len(tmp_list)

for gene in H:
    H.nodes[gene]['emp_p'] = sum(i>H.nodes[gene]['p_n'] for i in tmp_list)/float(total_scores)

H = permutation.multi_test_pvalues(H)

ennet.report(H, sys.argv[2])



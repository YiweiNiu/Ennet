#!/usr/bin/env python

import sys
from os import path

BINDIR = path.abspath(path.join(path.split(path.realpath(sys.argv[0]))[0], '..'))
sys.path.append(BINDIR)

import pickle
from Ennet import permutation
from Ennet import ennet

G = pickle.load(open(sys.argv[1], 'rb'))

H = G.copy()

tmp_dict = dict()
i = 0
for gene in H:
    for score in H.nodes[gene]['emp_p_n']:
        tmp_dict[i] = score
        i += 1

total_scores = i+1

for gene in H:
    H.nodes[gene]['emp_p'] = sum(tmp_dict[i]>H.nodes[gene]['p_n'] for i in tmp_dict)/float(total_scores)

H = permutation.multi_test_pvalues(H)

ennet.report(H, sys.argv[2])



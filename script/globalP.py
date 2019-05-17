#!/usr/bin/env python

import sys
import os
from os import path
import pickle

os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=1
os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=1
os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=1
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # export VECLIB_MAXIMUM_THREADS=1
os.environ["NUMEXPR_NUM_THREADS"] = "1" # export NUMEXPR_NUM_THREADS=1

import numpy as np

BINDIR = path.abspath(path.join(path.split(path.realpath(sys.argv[0]))[0], '..'))
sys.path.append(BINDIR)
from Ennet import permutation
from Ennet import ennet

G = pickle.load(open(sys.argv[1], 'rb'))

H = G.copy()

tmp_list = list()
for gene in H:
    tmp_list.extend(H.nodes[gene]['emp_p_n'])

tmp_list = np.array(tmp_list)
total_scores = len(tmp_list)

for gene in H:
    H.nodes[gene]['emp_p'] = np.sum(tmp_list>H.nodes[gene]['p_n'])/total_scores

H = permutation.multi_test_pvalues(H)

ennet.report(H, sys.argv[2])


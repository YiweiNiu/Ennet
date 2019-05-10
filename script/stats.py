#!/usr/bin/env python

import sys
from scipy.stats import hypergeom


def hypergeometric_test(A, B, population_size):
    """
    hypergeometri test

    @parameter A - a set of success of the population
    @parameter B - a set would be tested
    @parameter population_size - population size

    @return p-value of Hypergemetri Test
    """
    x = len(A & B)    # number of drawn "success"
    M = population_size    # population size
    n = len(A)    # number of success in the population
    N = len(B)    # sample size

    return hypergeom.sf(x-1, M, n, N)




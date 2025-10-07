
import numpy as np
import matplotlib.pyplot as plt
from numba import njit
import pandas as pd
import warnings
import pathlib




@njit  # just in time compilation speeds up sequence comparison, also see seq_compare_runtime file
def sc_jit(seq:bytes|str, ref:bytes|str):
    '''A fast Sequence Comparison of a short query sequence (seq) against a larger reference (ref). 
    Returns the number of matches for every position the sort sequence could bind to the larger one. Consideres only full overlap binding, i.e. no overhangs '''

    n = len(ref) - len(seq) + 1
    if n<1: print("sc_jit: wrong argument order")
    cv = np.zeros(n, dtype=np.uint16) 
    for i in range(n):
        matches = 0
        for j in range(len(seq)):
            if ref[i + j] == seq[j]:
                matches += 1
        cv[i] = matches
    return cv


@njit  # just in time compilation speeds up sequence comparison, also see seq_compare_runtime file
def sc_jit2(seq:bytes|str, ref:bytes|str):
    '''A fast Sequence Comparison of a short query sequence (seq) against a larger reference (ref).
    Returns the number of mismatches for every position the sort sequence could bind to the larger one. Consideres only full overlap binding, i.e. no overhangs '''

    n = len(ref) - len(seq) + 1
    if n<1: print("sc_jit: wrong argument order")
    cv = np.zeros(n, dtype=np.uint16) 
    for i in range(n):
        matches = 0
        for j in range(len(seq)):
            if ref[i + j] == seq[j]:
                matches += 1
        cv[i] = len(seq)-matches
    return cv

def random_seq(size:int):
    '''returns a random ACGT sequence of length = size'''
    return ''.join(np.array(["A", "C", "G", "T"])[np.random.randint(low = 0, high = 4, size = size)])

sc_jit(random_seq(10), random_seq(10));    # call it once so that numba compiles it




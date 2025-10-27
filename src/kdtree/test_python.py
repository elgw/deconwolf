#!/bin/env python

from sklearn.neighbors import NearestNeighbors
import numpy as np
import time
import sys


def print_memory():
    import os
    pid = os.getpid()
    pidfile = f"/proc/{pid}/status"
    with open(pidfile) as f:
        for line in f:
            if "VmPeak" in line:
                print(line.strip())

def test(n, k):
    rng = np.random.default_rng()
    X = rng.uniform([0, 0, 0], [1024, 1024, 1024], size=(n, 3))
    t0 = time.perf_counter()
    nbrs = NearestNeighbors(n_neighbors=k, algorithm='kd_tree').fit(X)
    t1 = time.perf_counter()
    # all vs all
    distances, indices = nbrs.kneighbors(X)
    t2 = time.perf_counter()
    print(f"Total {t2-t0:02f} s  Construction {t1-t0:02f} s  Query: {t2-t1:02f} s")

if __name__ == '__main__':
    n = 5000
    k = 5

    if( len(sys.argv) > 1):
        n = int(sys.argv[1]);
    if( len(sys.argv) > 2):
        k = int(sys.argv[2]);

    print(f"n = {n}, k = {k}")
    test(n, k)
    print_memory()

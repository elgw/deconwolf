#!/bin/env python

# To build the python wrapper for v 0.4.9
# cd src
# mkdir build
# cd build
# cmake ..
# make

progdesc = """Testing deconwolf via Python"""

import argparse
import sys
import pathlib

# if properly installed these two lines are not needed
here = pathlib.Path(__name__).parent.resolve()
sys.path.insert(1, str(here) + '/../src/build/')

import pydeconwolf as dw

print(f"deconwolf {dw.VERSION_MAJOR}.{dw.VERSION_MINOR}.{dw.VERSION_PATCH}")

import numpy as np
import matplotlib.pyplot as plt

def test_qalign(conf):

    nX = 90
    nY = 60
    X = np.random.random((nX, 3)).astype(np.float32)*1000
    Y = np.random.random((nY, 3)).astype(np.float32)*1000
    # This should fail
    try:
        H = dw.qalign(X=X, Y=Y)
    except:
        0

    nShared = 2 # We can find at most nS*(nS-1)/2 true correspondences
    D = 400*(np.random.random((1, 3))-0.5)
    X[0:nShared, :] = Y[0:nShared, :] + D
    shared = X[0:nShared, :].copy()
    np.random.shuffle(X)

    fig, ax = plt.subplots(ncols = 2, nrows = 2)
    ax[0][0].scatter(X[:, 0], X[:, 1], color='r', marker='x')
    ax[0][0].set_title("Point set A")
    ax[0][1].scatter(Y[:, 0], Y[:, 1], color='g', marker='o', facecolors='none')
    ax[0][1].set_title("Point set B")

    ax[1][0].scatter(X[:, 0], X[:, 1], color='r', marker='x')
    ax[1][0].scatter(Y[:, 0], Y[:, 1], color='g', marker='o', facecolors='none')
    for point in shared:
        ax[1][0].plot( [point[0], point[0]-D[0][0]], [point[1], point[1]-D[0][1]] , 'k')
    ax[1][0].set_title(f"A and B, the {nShared} shared points connected")
    if False:
        ax[0].scatter(shared[:, 0], shared[:, 1], color='b', marker='s', facecolors='none', s=100)
        ax[0].scatter(shared[:, 0] - D[0][0],
                  shared[:, 1] - D[0][1],
                  color='b', marker='s', facecolors='none', s=100)

    try:
        H = dw.qalign(X=X, Y=Y)
        delta = np.mean(H, axis=0)
        ax[1][1].scatter(X[:, 0], X[:, 1], color='r', marker='x')
        ax[1][1].scatter(Y[:, 0]+delta[0], Y[:, 1]+delta[1], color='g', marker='o', facecolors='none')
        ax[1][1].scatter(shared[:, 0], shared[:, 1], color='b', marker='s', facecolors='none', s=100)

        ax[1][1].set_title(f'qalign(A, b), delta={delta[0]:.1f}, {delta[1]:.1f}, {delta[2]:.1f}')
    except:
        ax[1][1].set_title(f'FAILED')

    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=progdesc,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--verbose', required=False, default=1, type=int,
                        help='Set the verbosity level')
    config = parser.parse_args()
    breakpoint()
    test_qalign(config)

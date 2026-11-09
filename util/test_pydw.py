# to build the python wrapper
# cd src
# mkdir build
# cd build
# cmake ..
# make



import sys
import pathlib

# if properly installed these two lines are not needed
here = pathlib.Path(__name__).parent.resolve()
sys.path.insert(1, str(here) + '/../src/build/')

import pydeconwolf as dw

print(f"deconwolf {dw.VERSION_MAJOR}.{dw.VERSION_MINOR}.{dw.VERSION_PATCH}")

import numpy as np
import matplotlib.pyplot as plt

def test_qalign():

    nX = 30
    nY = 60
    X = np.random.random((nX, 3)).astype(np.float32)*1000
    Y = np.random.random((nY, 3)).astype(np.float32)*1000
    # This should fail
    try:
        H = dw.qalign(X=X, Y=Y)
    except:
        0

    nShared = 3 # We can find at most nS*(nS-1)/2 true correspondences
    D = 200*(np.random.random((1, 3))-0.5)
    X[0:nShared, :] = Y[0:nShared, :] + D
    shared = X[0:nShared, :].copy()
    np.random.shuffle(X)

    fig, ax = plt.subplots(ncols = 2, nrows = 1)
    ax[0].scatter(X[:, 0], X[:, 1], color='r', marker='x')
    ax[0].scatter(Y[:, 0], Y[:, 1], color='g', marker='o', facecolors='none')
    ax[0].scatter(shared[:, 0], shared[:, 1], color='b', marker='s', facecolors='none', s=100)

    ax[0].scatter(shared[:, 0] - D[0][0],
                  shared[:, 1] - D[0][1],
                  color='b', marker='s', facecolors='none', s=100)

    ax[0].set_title(f'Input points, shared={nShared}')

    try:
        H = dw.qalign(X=X, Y=Y)
        delta = np.mean(H, axis=0)
        ax[1].scatter(X[:, 0], X[:, 1], color='r', marker='x')
        ax[1].scatter(Y[:, 0]+delta[0], Y[:, 1]+delta[1], color='g', marker='o', facecolors='none')
        ax[1].scatter(shared[:, 0], shared[:, 1], color='b', marker='s', facecolors='none', s=100)

        ax[1].set_title(f'Registered points delta={delta[0]:.1f}, {delta[1]:.1f}, {delta[2]:.1f}')
    except:
        ax[1].set_title(f'FAILED')

    plt.show()

test_qalign()

"""This script demonstrates how dw align-dots works.

It will try to align the points in the files A.dots.tsv to A.dots.tsv
using deconwolf and then plot the results.

To run the script set up with something like:

$ python -m venv .venv
$ source .venv/$ source .venv/bin/activate
$ pip install pandas matplotlib


You can also just run

$ dw align-dots a.dots.tsv b.dots.tsv

and see what happens.

Expected time consumption: < 0.1 s

Expected result:
156 correspondences were found after alignment
Shift: [-0.97, 75.04 -2.91] with high confidence, goodness = 1.00e+00

You can also decrease the number of poits that deconwolf uses, it
should work down to --npoint 17.

"""

import subprocess
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


cmd = ['dw', 'align-dots', 'A.dots.tsv', 'B.dots.tsv', '--out', 'results.csv']
print(f"Running {' '.join(cmd)}")
subprocess.run(cmd)

print("Reading the found displacement")
# see dw align-dots --help for info about the columns
R = pd.read_table('results.csv', header=None, delimiter=',')
D = np.array([R[0], R[1], R[2]])

print("Loading dots")
A = pd.read_table('A.dots.tsv')
B = pd.read_table('B.dots.tsv')
A = A.loc[A['value'] > 1000, :]
B = B.loc[B['value'] > 1000, :]


A2 = A.copy()
A2['f_x'] = A['f_x'] + D[0]
A2['f_y'] = A['f_y'] + D[1]
A2['f_z'] = A['f_z'] + D[2]


fig, ax = plt.subplots(nrows=1, ncols=2)
ax[0].scatter(A['f_x'], A['f_y'], marker='x')
ax[0].scatter(B['f_x'], B['f_y'], marker='.')
ax[0].legend({'A', 'B'})
ax[0].set_title('input')

ax[1].scatter(A2['f_x'], A2['f_y'], marker='x')
ax[1].scatter(B['f_x'], B['f_y'], marker='.')
ax[1].legend({'A2', 'B'})
ax[1].set_title('shifted, A2 = A + delta')

plt.show()


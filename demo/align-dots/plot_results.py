# Get started with something like:
# python -m venv .venv
# source .venv/$ source .venv/bin/activate
# pip install pandas matplotlib

import matplotlib.pyplot as plt
import pandas as pd

A = pd.read_table('A.dots.tsv')
B = pd.read_table('B.dots.tsv')
A = A.loc[A['value'] > 1000, :]
B = B.loc[B['value'] > 1000, :]


plt.scatter(A['f_x'], A['f_y'], marker='x')
plt.scatter(B['f_x'], B['f_y'], marker='.')
plt.legend({'A', 'B'})
plt.title('unshifted')
plt.show()


# run
# $ dw align-dots a.dots.tsv b.dots.tsv --radius 120
D = [-0.96,  75.06, -2.93]
A2 = A.copy()

A2['f_x'] = A['f_x'] + D[0]
A2['f_y'] = A['f_y'] + D[1]
A2['f_z'] = A['f_z'] + D[2]

plt.scatter(A2['f_x'], A2['f_y'], marker='x')
plt.scatter(B['f_x'], B['f_y'], marker='.')
plt.legend({'A2', 'B'})
plt.title('shifted, A2 = A + delta')
plt.show()


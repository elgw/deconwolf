#pip install scikit-learn
import numpy as np
from sklearn.datasets import *
import os

outdir = 'datasets/'

# https://stackoverflow.com/questions/43640546/how-to-make-randomforestclassifier-faster
datasets = ['iris', 'digits', 'wine', 'breast_cancer','diabetes', 'rand5_100', 'rand6_10']

if not os.path.exists(outdir):
    os.mkdir(outdir)


for name in datasets:
    outname = os.path.join(outdir, name + ".tsv")
    print(f"{outname}")

    if os.path.exists(outname):
        continue

    if name == 'rand5_100':
        ns = 100000
        nf = 100
        header = []
        for kk in range(0, nf):
            header.append(f'f{kk}')
        header.append('class')
        data = np.random.rand(ns, nf)
        target = np.random.randint(low=0, high=2, size=(ns, 1))
        data_array = np.concatenate((data, target), axis=1)
    elif name == 'rand6_10':
        ns = 1000000
        nf = 10
        header = []
        for kk in range(0, nf):
            header.append(f'f{kk}')
        header.append('class')
        data = np.random.rand(ns, nf)
        target = np.random.randint(low=0, high=2, size=(ns, 1))
        data_array = np.concatenate((data, target), axis=1)
    else:
        dataset = globals()['load_' + name]()
        data = dataset.data
        target = dataset.target
        target = target.reshape((len(target), 1))
        data_array = np.concatenate((data, target), axis=1)
        # breakpoint()
        header = dataset.feature_names
        header = list(map(str, header))
        header.append('class');


    np.savetxt(outname, data_array,
               header = '\t'.join(header),
               delimiter='\t', fmt='%.2f')

import time
import numpy as np
import os
import sys


def get_peak_memory():
    hwm = 0
    with open('/proc/self/status', 'r') as fid:
        for line in fid:
            if 'VmHWM' in line:
                parts = line.split(' ')
                hwm = int(parts[-2])

    return hwm



def read_tsv(fname):
    with open(fname, "r") as fid:
        header = fid.readline()
        data = np.fromfile(fid, sep='\t')
    header = header.split('\t')
    #breakpoint()
    data = data.reshape((len(data)//len(header), len(header)))
    return header, data


def print_peak_memory():
    with open('/proc/self/status', 'r') as fid:
        for line in fid:
            if 'VmPeak' in line:
                print(line.strip())
            if 'VmHWM' in line:
                print(line.strip())


if __name__ == '__main__':
    # print_peak_memory()

    (header, data) = read_tsv(sys.argv[1])
    X = data[:,0:-1]
    # Y = data[:,-1:]
    Y = data[:,-1]

    n_tree=100
    if len(sys.argv) > 2:
        n_tree = int(sys.argv[2])


    from sklearn.ensemble import RandomForestClassifier

    clf = RandomForestClassifier(n_estimators=n_tree)


    clf.n_jobs=-1
    clf.min_samples_split=2

    mem0 = get_peak_memory()
    t1 = time.perf_counter()
    clf = clf.fit(X, Y)
    t2 = time.perf_counter()
    mem1 = get_peak_memory()


    print(f"deltaRSS: {mem1-mem0} kb")

    t3 = time.perf_counter()
    P = clf.predict(X)
    t4 = time.perf_counter()


    print(f"sk: {np.sum(P==Y)} / {len(P)} predicted correctly")


    print(clf.get_params())

    print(f"tTraining: {t2-t1:.6f}")
    print(f"tPrediction {t4-t3:.6f}")

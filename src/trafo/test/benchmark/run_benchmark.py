import pandas as pd
import subprocess
import re
import io
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

datasets = [
    "iris",
    "digits",
    "wine",
    "breast_cancer",
    "diabetes"]

if True:
    datasets.append("rand5_100")
    datasets.append("rand6_10")

datadir = 'datasets'

def run_trafo(tsvfile, ntree=1):

    modelfile = os.path.join(datadir, tsvfile + ".trafo")

    ret = subprocess.run(['../../build/trafo_cli', '--train', os.path.join(datadir, tsvfile), '--cout', modelfile, "--ntree", str(ntree), '--verbose', '2'], capture_output=True)

    if len(ret.stderr) > 0:
        print(f'{ret.args} returned with the following error:')
        print(f'{ret.stderr}')
        exit(1)

    if ret.returncode != 0:
        print(f'{ret.args} returned an error but no message')
        exit(1)

    info = io.StringIO(ret.stdout.decode())

    t_train = -1
    mem = -1
    for line in info:
        if "trafo: Forest training took" in line:
            tmp = line.split()
            t_train = float(tmp[-2])
        if "deltaRSS" in line:
            tmp = line.split(' ')
            mem = float(tmp[-1])

    assert(t_train != -1)
    assert(mem != -1)

    ret = subprocess.run(['../../build/trafo_cli', '--predict', os.path.join(datadir, tsvfile), '--model', modelfile, '--verbose', '2'], capture_output=True)

    if len(ret.stderr) > 0:
        print(f'{ret.args} returned with the following error:')
        print(f'{ret.stderr}')
        exit(1)

    if ret.returncode != 0:
        print(f'{ret.args} returned an error but no message')
        exit(1)

    info = io.StringIO(ret.stdout.decode())

    t_predict = -1

    for line in info:
        if "trafo: Prediction took" in line:
            tmp = line.split()
            t_predict = float(tmp[-2])

    assert(t_predict != -1)

    return mem, t_train, t_predict


def run_skl(tsvfile, ntree=1):
    """ A new process is required in order to measure
    the memory usage
    """

    if ntree==1:
        ret = subprocess.run(['python', 'skl_rf_tree.py', os.path.join(datadir, tsvfile)], capture_output=True)
    else:
        ret = subprocess.run(['python', 'skl_rf_forest.py', os.path.join(datadir, tsvfile), str(ntree)], capture_output=True)

    if len(ret.stderr) > 0:
        print(f'{ret.args} returned with the following error:')
        print(f'{ret.stderr}')
        exit(1)

    if ret.returncode != 0:
        print(f'{ret.args} returned an error but no message')
        exit(1)


    info = io.StringIO(ret.stdout.decode())

    mem = -1
    tTrain = -1
    tPred = -1

    for line in info:
        if "deltaRSS" in line:
            mem = line.split(' ')
            mem = float(mem[1])
        if "tTraining" in line:
            tmp= line.split(' ')
            tTrain = float(tmp[1])
        if "tPrediction" in line:
            tmp = line.split(' ')
            tPred = float(tmp[1])

    assert(mem != -1)
    assert(tTrain != -1)
    assert(tPred != -1)
    return mem, tTrain, tPred


#
#
#

niter = 25
ntree = 100

if len(sys.argv) > 1:
    ntree = int(sys.argv[1])

print(f"ntree={ntree} niter={niter}")

df = pd.DataFrame(columns=['dataset','method','t_train', 't_predict', 'mem_fit_kb'])


for tsv in datasets:
    print(f"Dataset: {tsv} ", end="")
    for iter in range(0, niter):
        print(".", end ="")
        sys.stdout.flush()
        m, ttrain, tpredict = run_skl(tsv + '.tsv', ntree=ntree)
        row = {'dataset': [tsv], 'method': ['skl'], 't_train': [ttrain], 't_predict': [tpredict], 'mem_fit_kb': [m]}
        row = pd.DataFrame.from_dict(row)
        if(len(df) == 0):
            df = row
        else:
            df = pd.concat([df, row],  ignore_index=True)

        m, ttrain, tpredict = run_trafo(tsv + '.tsv', ntree=ntree)
        row = {'dataset': [tsv], 'method': ['trafo'], 't_train': [ttrain], 't_predict': [tpredict], 'mem_fit_kb': [m]}
        row = pd.DataFrame.from_dict(row)
        df = pd.concat([df, row],  ignore_index=True)
    print("")

if not os.path.exists('results'):
    os.mkdir(results)

outname = f"results/benchmarks_ntree{ntree}.tsv"
print(f"Writing to {outname}")
df.to_csv(outname, sep="\t")

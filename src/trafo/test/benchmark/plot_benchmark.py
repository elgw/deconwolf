import pandas as pd
import subprocess
import re
import io
import matplotlib.pyplot as plt
import numpy as np
import sys

df = pd.read_csv(sys.argv[1], sep="\t")


df['id'] = df['dataset'] + ":" + df['method']

# Sort by id
if 0:
    df.boxplot(column=['t_train'], by='id')
    plt.show()

    df.boxplot(column=['mem_fit_kb'], by='id')
    plt.show()

    df.boxplot(column=['t_predict'], by='id')
    plt.show()


# todo save table -- plot with separate script

# Todo: add titles etc...

# Better to create an new table, can be formated directly as markdown etc

ids = sorted(list(set(df['id'])))
print('Method          , t_train, mem_fit_kb')
summary = []
for id in ids:
    sub = df[df['id'] == id]
    print(f"{id:18s} {np.mean(sub['t_train']):.3f} {np.mean(sub['mem_fit_kb']):.0f}")

    row = {'dataset': [ sub['dataset'].iloc(0)[0] ],
           'method': [ sub['method'].iloc(0)[0] ],
           't_train_avg': [np.mean(sub['t_train'])],
           't_train_std': [np.std(sub['t_train'])],
           't_predict_avg': [np.mean(sub['t_predict'])],
           't_predict_std': [np.std(sub['t_predict'])],
           'mem_fit_kb': [np.mean(sub['mem_fit_kb'])],
           }
    row = pd.DataFrame.from_dict(row)
    if len(summary) == 0:
        summary = row
    else:
        summary = pd.concat([summary, row], ignore_index=True)

summary['mem_fit_kb'] = summary['mem_fit_kb'].astype(int)
print(summary.to_markdown(index=False)) # requires tabulate

skl = summary[summary['method'] == 'skl']
trafo = summary[summary['method'] == 'trafo']

np.array(skl['t_train_avg'])/np.array(trafo['t_train_avg'])
np.array(skl['t_predict_avg'])/np.array(trafo['t_predict_avg'])
np.array(skl['mem_fit_kb'])/np.array(trafo['mem_fit_kb'])


if 0:
    fig, ax = plt.subplots(1)
    ax.scatter(skl['t_train_avg'], trafo['t_train_avg'])
    ax.set_xlabel('scikit-learn')
    ax.set_ylabel('trafo')
    ax.set_title('avg t_train [s]')
    tmax = np.max(summary['t_train_avg'])
    ax.set_xlim([0, tmax])
    ax.set_ylim([0, tmax])
    ax.grid()
    fig.show()


breakpoint()

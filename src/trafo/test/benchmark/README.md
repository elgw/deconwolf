Extract some datasets to test on with

```
$ python extract_test_datasets.py
```

To measure the performance, use

```
$ python run_benchmark.py 1 # Test a single tree
# Generates: results/benchmarks_ntree1.tsv
$ python run_benchmark.py 100 # Test a forest with 100 trees
# Generates: results/benchmarks_ntree100.tsv
```

And summarize the results with command like

```
python plot_benchmark.py results/benchmarks_ntree100.tsv
```

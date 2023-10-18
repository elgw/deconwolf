Implements random forests for binary classification of pixel data.

Splitting nodes based on Gini's impurity.

## Notes to the next developer
 - Most time is spent sorting the features with qsort, for some data
   bucket sorting, count-sorting, etc, could make sense.
 - For creating training data, it seems like a good idea to use
   a standard drawing program like gimp, and put the class information in
   a separate layer.
 - Both Gini's impurity and entropy are special cases of Tsallis entropy.

## TODO
  - [ ] Refactor and start over again.
  - [ ] Extend the test suite... Constant input, invalid input, ...
 - [x] Convert tree to table.
 - [ ] Convert table to code.
 - [x] Make a node final if no threshold can be found.
 - [x] Parallelise training and classification with PrfForest.
 - [x] Boosting/bagging to create trees. Use `n_features = sqrt(q)`
        for each tree. Use `1 - 1/e` samples from the training set
        for each tree.
 - [x] Avoid malloc/free for individual nodes. Batch allocate and use buffer memory.
 - [x] Only consider `max_features = sqrt(q)` when consider splitting a node.
 - [x] Compile qsort with the program, increases speed with about 20 percent.
 - [ ] Just-In-Time (JIT) compilation of Trees, check out the [GCC JIT](https://gcc.gnu.org/onlinedocs/jit/intro/tutorial01.html)

## Files:
 - **prf_forest** random forest implementation
 - **prf_tree** decision tree implementation
 - **prf_ut** some unit tests, compile with **make prf_ut**
 - **prf_util** utility functions
 - **qsort** quick sort

## Timings:

### Dataset 1, `N = 150001`, `M = 14`
A Single tree, `minSize = 1`

| Method | Training | Classification |
| --- | --- | --- |
| `PrfTree` | 0.45 s | 0.0085 s |
| MATLAB/`fitctree` | 0.69 s | 0.043 s |

A Forest, `nTres = 200`, `minSize = 200`
| Method | Training | Classification |
| --- | --- | --- |
| `PrfForest` | 13.6 s | 1.42 s |
| MATLAB/`TreeBagger` | 57.2 | 5.2 |
| Python/sklearn | TODO | TODO |


# References:
 - [scikit-learn.org](https://scikit-learn.org/stable/modules/ensemble.html#random-forest-parameters)
 - [stackoverflow.com](https://stackoverflow.com/questions/40847745/subsample-size-in-scikit-learn-randomforestclassifier)

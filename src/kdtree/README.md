# k-d tree (but only the 3D case)

Version 0.1.1 2025-05-14.

The [K-d tree](https://en.wikipedia.org/wiki/K-d_tree) is a fun data
structure, useful for finding k-nearest neighbours and neighbours
within some distance in point clouds. The benefits it provides
compared to brute force drops quickly with the number of dimensions as
you can read on the Wiki page.

This repo supports exactly what I need and nothing more, so the
functionality is quite minimal. I'd be happy if anyone else finds it
useful and can send me a bug report now and then, or even a pull
request :)

A few short notes:
- No dependencies.
- 46 kb when compiled as a static library.
- Builds with gcc, clang and musl-gcc under linux.
- Usual warnings applies, use with caution!

## Usage
Below are examples of the supported methods:

``` C
#include <kdtree.h>
...
/* X: N 3D points [3 x N] */
kdtree_t * T = kdtree_new(X, N, 20);

/* Find the k nearest neighbours to Q [3 x 1] */
size_t * knn = kdtree_query_knn(T, Q, k);

/* Find any point within a distance of radius to Q */
size_t * idx = kdtree_query_radius(T, Q, radius, &n);

/* Evaluate the point density under the point, using
 * an isotropic Gaussian controlled by sigma.        */
double v = kdtree_kde(T, Q, sigma);

/* When done */
kdtree_free(T);
```

see `kdtree.h` for the complete function signatures and some
documentation. Look in `kdtree_ut.c` for complete usage examples.

## Details
- Data partitioning using Hoare's scheme, typically used in
  quicksort and quickselect.

- For finding the k nearest neighbours the candidates are put in a
  priority queue, implemented by a binary heap.

- The memory layout of the nodes has a big impact on performance and
  memory usage. The code in this repo use the
  [Eytzinger](https://arxiv.org/abs/1509.05053) layout, which is the
  same as used in [binary
  heaps](https://en.wikipedia.org/wiki/Binary_heap). On the positive
  side this give a good memory locality and fast queries. The major
  downside is that we need to decide upfront how deep the tree should
  be, which means that the memory usage (for the tree, excluding the
  data points) will grow in steps of approximately 2 when the number
  of points passes some boundaries.

- There is no parallel code for the tree construction at the moment
  although that would be possible to do. `kdtree_query_radius` and
  `kdtree_kde` should be thread safe. `kdtree_query_knn` is not tread safe
  but `kdtree_query_knn_multi` can be used to query multiple points in
  parallel.

- Can use GSL (`gsl_stats_median`) to find the pivot or the provided
  quick select implementation.

## Performance hints

Only using one thread. The function `kdtree_query_knn` found here is
denoted "this" in the table.

For reference, there are also results from
`sklearn.neighbors.NearestNeighbors`, see `test_python.py` for the
test code. Sklearn is probably an interface to
[ckdtree](https://github.com/scipy/scipy/tree/main/scipy/spatial/ckdtree/src)
but that is just a hypothesis, not a fact. The comparison is not fair,
comparisons seldom are, since the Python code stores the full
result (an Nxk matrix) at the end, the code here does not. For
sklearn, the memory measurement includes the whole Python environment,
not just the algorithm and the associated data.

Finding the k=5 nearest neighbours for each point among
N=5,000.

| Software | Tree construction |  Query | Total time |  VmPeak |
| -------- | ----------------- | ------ | ---------- | ------- |
| this     |            0.8 ms | 4.6 s  |     5.4 ms |    4 MB |
| sklearn  |            1.7 ms | 9.9 ms |    11.5 ms | 1504 MB |

N=1,000,000, k = 5

| Software | Tree construction | Query | Total time |  VmPeak |
| -------- | ----------------- | ----- | ---------- | ------- |
| this     |             0.3 s | 2.0 s |      2.3 s |   87 MB |
| sklearn  |             0.8 s | 5.8 s |      6.5 s | 1695 MB |


N=1,000,000, **k = 100**

| Software | Tree construction | Query | Total time |  VmPeak |
| -------- | ----------------- | ----- | ---------- | ------- |
| this     |            0.3 s  |17.2 s |     17.5 s |   87 MB |
| sklearn  |            0.8 s  |35.8 s |     36.6 s | 4664 MB |

N = 100,000,000, k = 5

| Software | Tree construction | Query | Total time |   VmPeak |
| -------- | ----------------- | ----- | ---------- | -------- |
| this     |              40 s | 483 s |      524 s |  8878 MB |
| sklearn  |             187 s | 968 s |     1155 s | 20580 MB |


## Current validation steps

More tests should be written. Especially to cover corner cases. The
current test pack includes:

- Compile with zero warnings using `gcc -Wall -Wextra -pedantic
  -std=gnu11 -g3 -fanalyzer`.
- Passes the few tests in `kdtree_ut.c`, some of them are comparisons
  to brute force calculations.
- **valgrind** finds no issues when running `kdtree_ut`.

## Build/Install

Use the makefile to generate the test program.

To build and install the library, please use

``` shell
mkdir build
cd build
cmake ..
make
sudo make install
```

To include in your project with CMake, you could copy the files to a
subfolder, for example named `modules/`, and then add something like
this to your `CMakeLists.txt`:

``` CMake
  add_subdirectory("modules/kdtree/")
  target_include_directories(myAPP PUBLIC "modules/kdtree/include/")
  target_link_directories(myAPP PUBLIC "modules/kdtree/")
  target_link_libraries(myAPP kdtree)
```

## To Do

- [ ] For N dimensions. Remove `#define KDTREE_DIM 3` and make it a
      parameter. Write tests and do the small adjustments needed.
- [ ] Fix so that `kdtree_query_knn` is thread safe -- via a query
      object containing the per-thread data?
- [ ] More validation.

## Changelog

- 0.1.1, 2025-04-05, added `kdtree_kde_mean` for mean shift algorithms.

- 0.1.0 2024-08-09.

## References

To read:

- [K-d tree page on Wikipedia](https://en.wikipedia.org/wiki/K-d_tree)
- [Blog post: Color quantization, minimizing variance, and k-d trees](https://www.crisluengo.net/archives/932/)

Implementations:

- Python: [sklearn.neighbors.KDTree](https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.KDTree.html)
- Matlab: [knnsearch](https://se.mathworks.com/help/stats/knnsearch.html)

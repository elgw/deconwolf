#ifndef __qsort_h_
#define __qsort_h_
#ifndef __compar_d_fn_t // Might be in stdlib.h
typedef int (*__compar_d_fn_t) (const void *, const void *, void *);
#endif

void
_quicksort (void *const pbase, size_t total_elems, size_t size,
	    __compar_d_fn_t cmp, void *arg);
#endif

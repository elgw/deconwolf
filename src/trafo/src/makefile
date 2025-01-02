# test: with some debug flags
# run: optimization flags


CFLAGS=-Wall -Wextra -pedantic -std=gnu11 -fopenmp
LDFLAGS=-lm -fopenmp

DEBUG?=1
ifeq ($(DEBUG),1)
	CFLAGS+=-g3
else
	CFLAGS+=-O3 -march=native -DNDEBUG
	LDFLAGS+=-flto
endif

FANALYZER?=0
ifeq ($(FANALYZER),1)
	CFLAGS+=-fanalyzer
endif

SAN?=0
ifeq ($(SAN),1)
	CFLAGS+=-fsanitize=address
endif

ifeq ($(SAN),0)
ifeq ($(FANALYZER),0)
ifeq ($(DEBUG), 1)
	CFLAGS+=-U_FORTIFY_SOURCE -D_FORTIFY_SOURCE -O
endif
endif
endif

# For usage with gprof
# ./run
# gprof run gmon.out > analysis

PG?=0
ifeq ($(PG),1)
	CFLAGS+=-pg
endif

FILES=qsort.c prf_util.c prf_forest.c test.c prf_tree.c
HEADERS=$(wildcard *.h)

LTRAFO_SRC = trafo.c sortbox.c qsort.c trafo_util.c ftab.c gini.c entropy.c
TRAFO_SRC = $(LTRAFO_SRC) trafo_cli.c
trafo_cli: $(TRAFO_SRC)
	$(CC) $(CFLAGS) $(TRAFO_SRC) $(LDFLAGS) -o trafo_cli

VISIBILITY=-DHAVE___ATTRIBUTE__VISIBILITY_HIDDEN -fvisibility=hidden

libtrafo.so: $(LTRAFO_SRC)
	$(CC) --shared -fPIC $(CFLAGS) $(VISIBILITY) $(LTRAFO_SRC) $(LDFLAGS) -o libtrafo.so

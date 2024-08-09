CC=gcc

CFLAGS=-Wall -Wextra -pedantic -std=gnu11 -Iinclude/
LDFLAGS=-lm -lpthread

GSL?=0

ifeq ($(GSL),1)
CFLAGS+=`pkg-config gsl --cflags` -DGSL
LDFLAGS+=`pkg-config gsl --libs`
endif

DEBUG?=0

ifeq ($(DEBUG),1)
CFLAGS+=-g3
else
CFLAGS+=-O3 -DNDEBUG
LDFLAGS+=-flto
endif

FANALYZER?=0
ifeq ($(FANALYZER),1)
CFLAGS+=-fanalyzer
endif

SRC=src/kdtree.c src/pqheap.c src/quickselect.c

kdtree_ut: $(SRC) src/kdtree_ut.c makefile
	$(CC) $(CFLAGS) $(SRC) src/kdtree_ut.c $(LDFLAGS) -o kdtree_ut

libkdtree.a: $(SRCFILES) makefile
	$(CC) -c $(CFLAGS) $(SRC) $(LDFLAGS)
	ar rcs libkdtree.a kdtree.o


libkdtree.so: $(SRCFILES) makefile
	$(CC) $(CFLAGS) -fPIC -shared $(SRC) $(LDFLAGS) -o libkdtree.so

install: libkdtree.so include/kdtree.h
	# TODO COPY to DESIRED PATHS

CFLAGS=-Wall -Wextra
CFLAGS+=-pedantic -std=gnu11 -fopenmp -std=gnu11
CFLAGS+=-O3 -DNDEBUG

LTRAFO_SRC=src/trafo.c \
src/qsort.c \
src/sortbox.c \
src/trafo_util.c \
src/ftab.c \
src/gini.c \
src/entropy.c

LTRAFO_OBJ=trafo.o \
qsort.o \
sortbox.o \
trafo_util.o \
ftab.o \
gini.o \
entropy.o

VISIBILITY=-DHAVE___ATTRIBUTE__VISIBILITY_HIDDEN -fvisibility=hidden


# $(OBJDIR)/%.o : %.c | $(OBJDIR)

%.o: src/%.c
	$(CC) $(CFLAGS) $< -c


libtrafo.a: $(LTRAFO_OBJ)
	ar -rc libtrafo.a *.o

libtrafo.so:
	$(CC) --shared -fPIC $(VISIBILITY) $(CFLAGS) \
$(LTRAFO_SRC)  -o libtrafo.so


clean:
	rm -rf *.o
	rm -rf *.a
	rm -rf *.so

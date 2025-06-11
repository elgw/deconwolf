CFLAGS=-Wall -Wextra
CFLAGS+=-pedantic -std=gnu11 -fopenmp -D_DEFAULT_SOURCE
DEBUG?=0

CFLAGS+=-fvisibility=hidden

ifeq ($(DEBUG),1)
CFLAGS += -O0 -Wno-unknown-pragmas -g3 -fanalyzer
else
CFLAGS += -O3 -DNDEBUG
endif

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

# For visibility of static library
CFLAGS+=-Wl,--exclude-libs=ALL

# $(OBJDIR)/%.o : %.c | $(OBJDIR)

%.o: src/%.c
	$(CC) $(CFLAGS) $< -c


libtrafo.a: $(LTRAFO_OBJ)
	ar -rc libtrafo.a *.o

libtrafo.so:
	$(CC) --shared -fPIC $(VISIBILITY) $(CFLAGS) \
$(LTRAFO_SRC)  -o libtrafo.so

CLI_FILES = $(wildcard src/*.c)
LDFLAGS=-lm

trafo_cli: $(CLI_FILES)
	$(CC) $(CFLAGS) $(CLI_FILES) $(LDFLAGS) -o trafo_cli

clean:
	rm -rf *.o
	rm -rf *.a
	rm -rf *.so

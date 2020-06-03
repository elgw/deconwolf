# Normal build:
# make -B
# To build with debug flags:
# make DEBUG=1 -B
# To use openMP:
# make OMP=1 -B

CC_VERSION = "$(shell gcc --version | head -n 1)"
GIT_VERSION = "$(shell git log --pretty=format:'%aD:%H' -n 1)"

XFLAGS = -DCC_VERSION=\"$(CC_VERSION)\"
XFLAGS += -DGIT_VERSION=\"$(GIT_VERSION)\"

DEBUG?=0
ifeq ($(DEBUG),1)
    CFLAGS =-Wall -g3 -DDEBUG 
else
    CFLAGS=-Wall -Wno-unknown-pragmas -DNDEBUG -O3 -flto -march=native -ftree-vectorize
endif

OMP?=0
ifeq ($(OMP), 1)
  CFLAGS += -fopenmp
else
  CFLAGS += -fno-openmp
endif

CFLAGS += $(XFLAGS)

CC = cc $(CFLAGS)

EXECUTABLE = bin/deconwolf
OBJECTS = fim.o tiling.o fft.o fim_tiff.o dw.o deconwolf.o
LIBRARIES = -lm -lfftw3f -lfftw3f_threads -ltiff
SRCDIR = src/

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) -o $@ $^ $(LIBRARIES)

%.o: $(SRCDIR)%.c
	$(CC) -c $<

clean:
	rm -f $(EXECUTABLE) $(OBJECTS)

install:
	cp bin/deconwolf /usr/bin/
	cp doc/deconwolf.1 .
	gzip deconwolf.1
	mv deconwolf.1.gz /usr/share/man/man1/

uninstall:
	rm /usr/bin/deconwolf
	rm /usr/share/man/man1/deconwolf.1.gz

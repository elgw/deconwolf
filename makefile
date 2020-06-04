# Normal build:
# make -B
# To build with debug flags and no OpenMP
# make DEBUG=1 OMP=0 -B
# For normal build
# make -B

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

OMP?=1
ifeq ($(OMP), 1)
  CFLAGS += -fopenmp
else
  CFLAGS += -fno-openmp
endif

CFLAGS += $(XFLAGS)

CC = cc $(CFLAGS)
SRCDIR = src/

dw = bin/deconwolf 
dw_OBJECTS = fim.o tiling.o fft.o fim_tiff.o dw.o deconwolf.o
dw_LIBRARIES = -lm -lfftw3f -lfftw3f_threads -ltiff

dwtm = bin/dw_tiffmax
dwtm_OBJECTS = fim_tiff.o deconwolf_tif_max.o
dwtm_LIBRARIES = -lm -ltiff


all: $(dw) $(dwtm)

$(dwtm): $(dwtm_OBJECTS)
	$(CC) -o $@ $^ $(dwtm_LIBRARIES)

$(dw): $(dw_OBJECTS)
	$(CC) -o $@ $^ $(dw_LIBRARIES)

%.o: $(SRCDIR)%.c
	$(CC) -c $<

clean:
	rm -f $(dw) $(dw_OBJECTS)
	rm -f $(dwtm) $(dwtm_OBJECTS)

install:
	cp bin/deconwolf /usr/bin/
	cp bin/dw_tiffmax /usr/bin
	cp doc/deconwolf.1 .
	gzip deconwolf.1
	mv deconwolf.1.gz /usr/share/man/man1/

uninstall:
	rm /usr/bin/deconwolf
	rm /usr/bin/dw_tiffmax
	rm /usr/share/man/man1/deconwolf.1.gz


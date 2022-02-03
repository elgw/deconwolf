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

CFLAGS = -Wall -Wextra -std=gnu99 -march=native -mtune=native
# Since we don't check for these errors, see man math_error
# we can disable this feature
CFLAGS += -fno-math-errno

DESTDIR?=/usr/local/bin
DEBUG?=0
ifeq ($(DEBUG),1)
    CFLAGS += -g3 -DDEBUG
else
   #CFLAGS +=  -g -O2 -ftree-vectorize -Wno-unknown-pragmas -flto
   CFLAGS +=  -O3 -Wno-unknown-pragmas -flto -DNDEBUG
   #-fno-math-errno no relevant performance gain
endif


dw_LIBRARIES =  -lm -ltiff
dwtm_LIBRARIES =  -lm -ltiff
dwbw_LIBRARIES = -lm -ltiff -lpthread -ltiff  -lgsl -lgslcblas

MKL?=0
ifeq ($(MKL),1)
CFLAGS += -DMKL `pkg-config mkl-static-lp64-iomp --cflags`
dw_LIBRARIES += `pkg-config mkl-static-lp64-iomp --cflags --libs`
dwtm_LIBRARIES += `pkg-config mkl-static-lp64-iomp --cflags --libs`
dwbw_LIBRARIES += `pkg-config mkl-static-lp64-iomp --cflags --libs`
else
dw_LIBRARIES +=  -lfftw3f
dwtm_LIBRARIES += -ltiff -lfftw3f
dwbw_LIBRARIES += -lfftw3f
endif
MANPATH=/usr/share/man/man1/

# on MacOS add -Xpreprocessor
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    CFLAGS +=
ifeq ($(MKL),0)
    dw_LIBRARIES += -lfftw3f_threads
endif
endif

ifeq ($(UNAME_S),Darwin)
    CFLAGS += -Xpreprocessor -fopenmp
    dw_LIBRARIES += -lomp -lfftw3f_threads
    MANPATH=/usr/local/share/man/man1
endif

ifeq ($(WINDOWS),1)
	CFLAGS += -DWINDOWS
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

dw = bin/dw
dw_OBJECTS = fim.o tiling.o fft.o fim_tiff.o dw.o deconwolf.o deconwolf_tif_max.o

dwbw = bin/dw_bw
dwbw_OBJECTS = fim.o fim_tiff.o dw_bwpsf.o bw_gsl.o lanczos.o li.o

# dwtm = bin/dw_tiffmax
# dwtm_OBJECTS = fim.o fim_tiff.o deconwolf_tif_max.o

all: $(dw) $(dwtm) $(dwbw)
# all: $(dw) $(dwtm) $(dwbw)


$(dw): $(dw_OBJECTS)
	$(CC) -o $@ $^ $(dw_LIBRARIES)

$(dwbw): $(dwbw_OBJECTS)
	$(CC) -o $@ $^ $(dwbw_LIBRARIES)

%.o: $(SRCDIR)%.c
	$(CC) -c $<

clean:
	rm -f $(dw) $(dw_OBJECTS)
	rm -f $(dwbw) $(dwbw_OBJECTS)

install:
	# Binaries
	cp bin/dw_bw $(DESTDIR)/dw_bw
	cp bin/dw $(DESTDIR)/dw
	# Man pages
	cp doc/dw.1 $(MANPATH)/dw.1
	cp doc/dw.1 $(MANPATH)/deconwolf.1
	cp doc/dw_bw.1 $(MANPATH)/dw_bw.1


uninstall:
	rm $(DESTDIR)/dw
	rm $(DESTDIR)/dw_bw
	rm $(DESTDIR)/dw_batch
	rm $(MANPATH)/dw.1
	rm $(MANPATH)/deconwolf.1
	rm $(MANPATH)/dw_bw.1
